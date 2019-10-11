#!/usr/bin/python
import gpxpy
import sys
from math import sin
from math import cos
from math import atan2
from math import sqrt
from math import radians
import argparse

parser = argparse.ArgumentParser(description='Calculate gain/loss segments for PCT.')
parser.add_argument('--segment_size',
        type=float,
        default=1.0,
        help='Size of individual segment in miles.')
parser.add_argument('--true_miles',
        type=float,
        default=0.0,
        help=('The known true distance of the route. If specified, raw gpx distances'
            ' will be adjusted to match the known distance. The scaling is linear.'))
parser.add_argument('--true_gain',
	type=float,
	default=0.0,
	help=('The known true elevation gain of the segment. Calculated gain/loss'
            ' will be adjusted to this number.'))
parser.add_argument('--true_loss',
        type=float,
        default=0.0,
        help=('The known true elevation loss of the segment. Calculated gain/loss'
            ' will be adjusted to this number.'))
parser.add_argument('--smoothing',
        type=float,
        default=0.3,
        help=('When calculating gain/loss, smooth neighboring gpx points by this'
            ' factor. This helps reduce the noisiness of the data. Let S be the'
            ' smoothing factor. Then the elevation of gpx point P(n) is'
            ' recalculated as follows:'
            ' P_Smoothed(n) := (1 - S) * P(n) + S / * (P(n-1) + P(n+1))'))
parser.add_argument(
        'gpxfile',
        metavar='N',
        type=str,
        nargs='+',
        help='gpx files to be processed')
parser.add_argument(
        '--section_header',
        type=str,
        default='',
        help='If set, print out section headers.')
# TODO(jaro): note that processing multiple gpx files don't necessarily make sense with --true_miles

args = parser.parse_args()


# Diameter in meters
R_EARTH_METERS = 6367 * 1000
R_EARTH_MILES = 3956

# Kilometers to miles conversion
# 6367 km = 3956 mi ---> 6367 / 3956 = km / mi
MI_PER_KM = 3956.0 / 6367.0
FEET_PER_METER = 3.28084

def lldist_haversine(x, y):
	"""Calculates distance between two gpx points in miles."""
	dlon = radians(y.longitude) - radians(x.longitude)
	dlat = radians(y.latitude) - radians(x.latitude)
	a = sin(dlat/2)**2 + cos(radians(x.latitude)) * cos(radians(y.latitude)) * (sin(dlon/2)**2)
	c = 2 * atan2(sqrt(a), sqrt(1-a))
	d = R_EARTH_MILES * c
	return d

def to_cart(p, alt = R_EARTH_MILES):
	"""Converts from polar to cartesian, using elevation in miles."""
	x = alt * cos(radians(p.latitude)) * sin(radians(p.longitude))
	y = alt * sin(radians(p.latitude))
	z = alt * cos(radians(p.latitude)) * cos(radians(p.longitude))
	return (x,y,z)

def lldist_cartesian(p1,p2):
	"""Calculates distance using polar->cartesian conversion and euclidean dist."""
	c1 = to_cart(p1)
	c2 = to_cart(p2)
	dist = 0
	for a,b in zip(c1,c2):
		dist += (b-a) ** 2
	return sqrt(dist)

def lldist_cartesian_elev(p1,p2):
	# Meters to miles = meters / 1000 / MI_PER_KM
	c1 = to_cart(p1, alt=R_EARTH_MILES + p1.elevation / 1000.0 / MI_PER_KM)
	c2 = to_cart(p2, alt=R_EARTH_MILES + p2.elevation / 1000.0 / MI_PER_KM)
	dist = 0
	for a,b in zip(c1,c2):
		dist += (b-a) ** 2
	return sqrt(dist)



def distance(points, dist_fn=lldist_haversine):
    """Calculates total distance between points in miles."""
    return sum(dist_fn(p,q) for p,q in zip(points, points[1:]))


# Calculating elevation gain can be easy/complex based on the desired precision.
# GPS data is inherently unstable/noisy. We can simply interpolate the three neighboring
# points to smooth out the variability and then calculate deltas among those.
def elevation_gain_loss(points):
    """Returns (gain, loss) tuple (in meters). Smooths out elevations of neighboring pts.

    smoothing: Number between 0.0 and 1.0 determining the weight of the point itself when
      averaging elevation measurements. Remainder of 1-smoothing will be assigned to the
      points previous/next neighbors distributed equally among the two.
    """
    smoothing = args.smoothing
    elevations = [p.elevation for p in points]
    smoothed = []
    for i, x in enumerate(elevations):
        prev = x
        next = x
        if i > 0:
            prev = elevations[i-1]
        if i < len(elevations) - 1:
            next = elevations[i+1]
        smoothed.append((1 - smoothing) * x + 0.5 * smoothing * (prev + next))

    #for x in elevations:
    #    avg_over = smoothed[-2:] + [x]  # Last 2 datapoints + this one
    #    smoothed.append(float(sum(avg_ovrr)) / len(avg_over))
    deltas = [b-a for a,b in zip(smoothed, smoothed[1:])]
    return (sum(filter(lambda x: x>0, deltas)), sum(filter(lambda x: x<0, deltas)))

# Pick distance measurement.
lldist = lldist_cartesian

total_dist_h = 0
total_dist_c = 0
total_dist_ce = 0
prev_point = None
total_gain = 0
total_loss = 0
points = 0
gain, loss = 0, 0
m2ft = lambda x: x * FEET_PER_METER


def get_pct_segment(filename):
    gpx = gpxpy.parse(open(filename, "r"))
    return gpx.tracks[0].segments[0]


# The goal is now to turn the series of gps points into series of segments with
# fixed distances and calculating (pct_mile_marker_start, total_gain_ft, total_loss_ft)
# Granularity of the output should be controlled by segment_size (miles). Approximation
# when two points are not exactly matching the boundaries may be necessary (or not).

def smoothed_elevation_ft(idx, points, smoothing=0.3):
    """Calculates smoothed elevation of idx-th point in points."""
    curr = points[idx]
    prev = points[idx]
    next = points[idx]
    if idx > 0:
        prev = points[idx-1]
    if idx < len(points) - 1:
        next = points[idx+1]

    return m2ft((1 - smoothing) * curr.elevation
                + 0.5 * smoothing * (prev.elevation + next.elevation))


def make_segments(points, segment_size=1, dist_fn=lldist_haversine, true_miles=None):

    """Returns list of (segment_end_miles, segment_gain_ft, segment_loss_ft, num_pts).

    Args:
      true_miles: if set, this specifies the actual path distance in miles.
        this will be used to even out the actual distance measurements to add up to
        the expected distance.
    """
    distance_mult = 1.0
    if true_miles:
        measured_miles = distance(points)
        # The points should be multiplied by true_miles / measured_miles to end up
        # with true_miles at the end.
        distance_mult = float(true_miles) / measured_miles
    total_miles = 0
    total_gain = 0
    total_loss = 0
    segment_start_miles = 0
    segment_gain = 0
    segment_loss = 0
    prev_ft = None
    prev_point = None
    num_points = 0
    for i, point in enumerate(points):
        # Initialize the previous point and move on
        if prev_point is None:
            prev_point = point
            prev_ft = smoothed_elevation_ft(i, points)
            continue

        num_points += 1
        total_miles += distance_mult * dist_fn(prev_point, point)
        curr_ft = smoothed_elevation_ft(i, points)
        if curr_ft > prev_ft:
            total_gain += curr_ft - prev_ft
            segment_gain += curr_ft - prev_ft
        else:
            total_loss += curr_ft - prev_ft
            segment_loss += curr_ft - prev_ft

        # Move the needle
        prev_point = point
        prev_ft = curr_ft

        # Segment is "at least" segment_size now, emit
        if total_miles - segment_start_miles >= segment_size:
            yield (total_miles, segment_gain, segment_loss, num_points)
            segment_start_miles = total_miles
            segment_gain = 0
            segment_loss = 0
            num_points = 0

    if total_miles > segment_start_miles:
        yield (total_miles, segment_gain, segment_loss, num_points)



# print "Section\tMiles\tGain(ft)\tLoss(ft)"
for f in args.gpxfile:
    s = get_pct_segment(f)
    #miles = distance(s.points)
    #gain, loss = elevation_gain_loss(s.points)
    #print "%s\t%.2f\t%.2f\t%.2f" % (f, miles, m2ft(gain), m2ft(loss))
    # For spreadsheet processing, get rid of unnecessary values
    # print "%.2f\t%.2f" % (m2ft(gain), m2ft(loss))
    # Now print the segmentations
    # print "Total points in segment: {0}".format(len(s.points))
    sg = make_segments(s.points, segment_size=args.segment_size, true_miles=args.true_miles)
    sg_gain = 0
    sg_loss = 0
    header = ''
    if args.section_header:
        header = '%s\t' % args.section_header
    for s in sg:
        print (header + '%2f\t%.2f\t%.2f\t%d' % s)
        sg_gain += s[1]
        sg_loss += s[2]
    # print '*** Total gain %.2f loss %.2f' % (sg_gain, sg_loss)
