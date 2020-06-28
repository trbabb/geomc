#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Shape

// #include <iostream>

#include <boost/test/unit_test.hpp>
#include <geomc/shape/Oriented.h>
#include <geomc/shape/Cylinder.h>


using namespace geom;
using namespace std;



BOOST_AUTO_TEST_SUITE(shape)


BOOST_AUTO_TEST_CASE(create_oriented) {
    typedef Oriented<Cylinder<double,3>> ocyl_t;
    ocyl_t ocyl = ocyl_t(Cylinder<double,3>());
    BOOST_CHECK(ocyl.contains(Vec3d()));
    ocyl.convexSupport(Vec3d(1,0,0));
}


BOOST_AUTO_TEST_SUITE_END()
