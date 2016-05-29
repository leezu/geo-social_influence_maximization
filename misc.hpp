#pragma once

namespace gsinfmax {
	namespace misc {
		double to_radians(double degree);
		double great_circle_length(double lat1, double lon1, double lat2, double lon2);
	}
}
