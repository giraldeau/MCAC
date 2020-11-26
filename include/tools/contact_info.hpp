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
    std::weak_ptr<Sphere> other_sphere;
    HalfSphereListContactInfo() = default;
    explicit HalfSphereListContactInfo(const std::shared_ptr<Sphere>& _other_sphere) :
        other_sphere(_other_sphere) {
    }
    HalfSphereListContactInfo(SphereContactInfo contact_info, const std::shared_ptr<Sphere>& _other_sphere) :
        SphereContactInfo(contact_info),
        other_sphere(_other_sphere) {
    }
};

class SphereListContactInfo : public HalfSphereListContactInfo {
public:
    std::weak_ptr<Sphere> moving_sphere;
    SphereListContactInfo() = default;
    SphereListContactInfo(HalfSphereListContactInfo contact_info, const std::shared_ptr<Sphere>&  _moving_sphere) :
        HalfSphereListContactInfo(contact_info),
        moving_sphere(_moving_sphere) {
    }
    SphereListContactInfo(SphereContactInfo contact_info, const std::shared_ptr<Sphere>& _moving_sphere, const std::shared_ptr<Sphere>& _other_sphere) :
        HalfSphereListContactInfo(contact_info, _other_sphere),
        moving_sphere(_moving_sphere) {
    }
    SphereListContactInfo(const std::shared_ptr<Sphere>& _moving_sphere, const std::shared_ptr<Sphere>& _other_sphere) :
        HalfSphereListContactInfo(_other_sphere),
        moving_sphere(_moving_sphere) {
    }
};

class AggregateContactInfo : public SphereListContactInfo {
public:
    std::weak_ptr<Aggregate> moving_aggregate;
    std::weak_ptr<Aggregate> other_aggregate;
    AggregateContactInfo() = default;
    AggregateContactInfo(SphereListContactInfo contact_info,
                         const std::shared_ptr<Aggregate>& _moving_aggregate,
                         const std::shared_ptr<Aggregate>& _other_aggregate) :
        SphereListContactInfo(contact_info),
        moving_aggregate(_moving_aggregate),
        other_aggregate(_other_aggregate) {
    }
};
}  // namespace mcac


#endif //INCLUDE_TOOLS_CONTACT_INFO_HPP


