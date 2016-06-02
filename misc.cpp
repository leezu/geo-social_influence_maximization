#include "misc.hpp"

#include <cmath>

namespace gsinfmax {
namespace misc {
double to_radians(double degree) { return degree * M_PI / 180; }

double great_circle_length(double lat1, double lon1, double lat2, double lon2) {
    double r{6378137.0}; // earth radius in meter

    lat1 = to_radians(lat1);
    lon1 = to_radians(lon1);
    lat2 = to_radians(lat2);
    lon2 = to_radians(lon2);

    double dlat = lat2 - lat1;
    double dlon = lon2 - lon1;

    return 2 * r * asin(sqrt(pow(sin(dlat / 2), 2) +
                             cos(lat1) * cos(lat2) * pow(sin(dlon / 2), 2)));
}
}
}
