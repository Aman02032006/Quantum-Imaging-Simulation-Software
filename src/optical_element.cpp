#include "optical_element.hpp"

// Constructor
OpticalElement::OpticalElement(const vec3 &pos, const vec3 &orient, const std::string &n)
    : position(pos), orientation(unit_vector(orient)), name(n) {}

// Getters
vec3 OpticalElement::getPosition() const
{
    return position;
}

vec3 OpticalElement::getOrientation() const
{
    return orientation;
}

std::string OpticalElement::getName() const
{
    return name;
}
