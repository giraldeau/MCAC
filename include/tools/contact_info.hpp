/*
 * MCAC
 * Copyright (C) 2020 CORIA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INCLUDE_TOOLS_CONTACT_INFO_HPP
#define INCLUDE_TOOLS_CONTACT_INFO_HPP
namespace mcac {
class SphereContactInfo {
public:
    double distance;
    SphereContactInfo() :
        distance(std::numeric_limits<double>::infinity()) {
    }
    explicit SphereContactInfo(double _distance) :
        distance(_distance) {
    }
    bool operator<=(const double _distance) const {
        return distance <= _distance;
    }
    bool operator<(const SphereContactInfo &b) const {
        return distance < b.distance;
    }
    explicit operator double() const {
        return distance;
    }
};

class HalfSphereListContactInfo : public SphereContactInfo {
public:
    size_t other_sphere{0};
    HalfSphereListContactInfo() = default;
    explicit HalfSphereListContactInfo(size_t _other_sphere) :
        other_sphere(_other_sphere) {
    }
    HalfSphereListContactInfo(SphereContactInfo contact_info, size_t _other_sphere) :
        SphereContactInfo(contact_info),
        other_sphere(_other_sphere) {
    }
};

class SphereListContactInfo : public HalfSphereListContactInfo {
public:
    size_t moving_sphere{0};
    SphereListContactInfo() = default;
    SphereListContactInfo(HalfSphereListContactInfo contact_info, size_t _moving_sphere) :
        HalfSphereListContactInfo(contact_info),
        moving_sphere(_moving_sphere) {
    }
    SphereListContactInfo(SphereContactInfo contact_info, size_t _moving_sphere, size_t _other_sphere) :
        HalfSphereListContactInfo(contact_info, _other_sphere),
        moving_sphere(_moving_sphere) {
    }
    SphereListContactInfo(size_t _moving_sphere, size_t _other_sphere) :
        HalfSphereListContactInfo(_other_sphere),
        moving_sphere(_moving_sphere) {
    }
};

class AggregateContactInfo : public SphereListContactInfo {
public:
    size_t moving_aggregate{0};
    size_t other_aggregate{0};
    AggregateContactInfo() = default;
    AggregateContactInfo(SphereListContactInfo contact_info, size_t _moving_aggregate, size_t _other_aggregate) :
        SphereListContactInfo(contact_info),
        moving_aggregate(_moving_aggregate),
        other_aggregate(_other_aggregate) {
    }
};
}  // namespace mcac


#endif //INCLUDE_TOOLS_CONTACT_INFO_HPP


