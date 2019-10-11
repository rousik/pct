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
        default=0.0,
        help=('When calculating gain/loss, smooth neighboring gpx points by this'
            ' factor. This helps reduce the noisiness of the data. Let S be the'
            ' smoothing factor. Then the elevation of gpx point P(n) is'
            ' recalculated as follows:'
            ' P_Smoothed(n) := (1 - S) * P(n) + S / * (P(n-1) + P(n+1))'))
parser.add_argument('--gain_threshold',
        type=int,
        default=5,
        help=('When calculating gain/loss, only consider points if their relative'
            ' difference (in meters) exceeds this threshold.'))
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


def smooth_elevations(elevations, smoothing=.25):
    smoothed = []
    for i, x in enumerate(elevations):
        prev = x
        next = x
        if i > 0:
            prev = elevations[i-1]
        if i < len(elevations) - 1:
            next = elevations[i+1]
        smoothed.append((1 - smoothing) * x + 0.5 * smoothing * (prev + next))
    return smoothed

# Calculating elevation gain can be easy/complex based on the desired precision.
# GPS data is inherently unstable/noisy. We can simply interpolate the three neighboring
# points to smooth out the variability and then calculate deltas among those.
def elevation_gain_loss(points):
    """Returns (gain, loss) tuple (in meters). Smooths out elevations of neighboring pts.
    """
    elevations = [p.elevation for p in points]
    if args.smoothing:
        elevations = smooth_elevations(elevations, smoothing=args.smoothing)

    gain = 0
    loss = 0
    last_e = elevations[0]
    for e in elevations:
        diff = e - last_e
        if abs(diff) < args.gain_threshold:
            continue
        last_e = e
        if diff > 0:
            gain += diff
        else:
            loss += diff
    return (gain, loss)

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
def split_by_distance(points, segment_length=1, dist_fn=lldist_haversine, true_miles=0.0):
    """This function splits list of points into list of sub-list where each sublist has specified distance.

    Args:
        segment_length: expected distance of each sub-list/segment in miles.
        dist_fn: function to use when calculating distance between gpx points.
        true_miles: Known distance of the entire segment. In case exact distance
          measurements are off by some number, the overall distances will be
          corrected so that the total distance matches true_miles.
    """
    distance_mult = 1.0
    if true_miles:
        measured_miles = distance(points)
        # The points should be multiplied by true_miles / measured_miles to end up
        # with true_miles at the end.
        distance_mult = float(true_miles) / measured_miles

    output_segments = []
    curr_segment = [] 
    curr_segment_length = 0
    last_pt = points[0]
    for pt in points:
        dist = distance_mult * dist_fn(last_pt, pt)
        # TODO(jaro): we could improve accurracy by:
        # 1. picking to include this point into previous/next segment based on
        #    the distance to the boundary point (this may not make a difference in
        #    aggregate)
        # 2. constructing split-points and averaging their elevations.
        #
        # As of now, we simply add up points until they exceed the segment length
        # and average these numbers afterwards.
        curr_segment_length += dist
        curr_segment.append(pt)
        last_pt = pt

        if curr_segment_length >= segment_length:
            yield curr_segment
            curr_segment = [pt]
            curr_segment_length = 0

    # Add the last (possibly shorter) segment
    if curr_segment_length > 0:
        yield curr_segment


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
    dist_mult = 1.0
    if args.true_miles:
        dist_mult = args.true_miles / distance(s.points)

    segs = split_by_distance(s.points, segment_length=args.segment_size, true_miles=args.true_miles)

    header = ''
    if args.section_header:
        header = '%s\t' % args.section_header

    total_gain, total_loss = 0, 0
    total_distance = 0
    for sg in segs:
        sg_gain, sg_loss = elevation_gain_loss(sg)
        total_gain += sg_gain
        total_loss += sg_loss
        total_distance += dist_mult * distance(sg)
        print (header + '%2f\t%.2f\t%.2f\t%d' % (
            total_distance, m2ft(sg_gain), m2ft(sg_loss), len(sg)))
    print '*** Total gain for %s %.2f loss %.2f' % (
            f,
            m2ft(total_gain), m2ft(total_loss))
