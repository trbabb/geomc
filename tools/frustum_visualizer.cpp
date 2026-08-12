#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <geomc/shape/Frustum.h>
#include <geomc/shape/Rect.h>
#include <geomc/shape/Simplex.h>
#include <geomc/shape/Sphere.h>

using namespace geom;

namespace {

using T = double;
using P2 = Vec<T,2>;
using P3 = Vec<T,3>;
using Segment = std::pair<P3,P3>;

struct Options {
    std::string output = "frustum_samples.json";
    std::vector<std::string> shapes {"rectangle", "circle", "triangle"};
    std::uint64_t seed = std::random_device{}();
    int samples = 64;
    int rays = 1000;
};

struct Sample {
    std::string shape;
    std::string description;
    Rect<T,1> height;
    P3 center;
    T radius;
    std::vector<Segment> wire;
    std::vector<Segment> outside;
    std::vector<Segment> hits;
};

T uniform(std::mt19937_64& rng, T lo, T hi) {
    return std::uniform_real_distribution<T>(lo, hi)(rng);
}

P3 random_unit(std::mt19937_64& rng) {
    std::normal_distribution<T> normal;
    P3 p;
    do {
        p = {normal(rng), normal(rng), normal(rng)};
    } while (p.mag2() == 0);
    return p.unit();
}

Rect<T,1> random_height(std::mt19937_64& rng) {
    T a = uniform(rng, 0.15, 1.2);
    T b = uniform(rng, 1.3, 3.5);
    switch (std::uniform_int_distribution<int>(0, 2)(rng)) {
        case 0: return {a, b};
        case 1: return {-b, -a};
        default: return {-a, b}; // Frustum clips this to [0,b].
    }
}

void add_wireframe(std::vector<Segment>& wire, const std::vector<P2>& boundary,
                   Rect<T,1> height) {
    for (T h : {height.lo, height.hi}) {
        for (std::size_t i = 0; i < boundary.size(); ++i) {
            P2 a = h * boundary[i];
            P2 b = h * boundary[(i + 1) % boundary.size()];
            wire.push_back({P3{a, h}, P3{b, h}});
        }
    }
    for (const P2& p : boundary) {
        P2 a = height.lo * p;
        P2 b = height.hi * p;
        wire.push_back({P3{a, height.lo}, P3{b, height.hi}});
    }
}

template <typename Base>
Sample make_sample(const std::string& name, const std::string& description,
                   const Base& base, const std::vector<P2>& boundary,
                   Rect<T,1> requested_height, std::mt19937_64& rng, int ray_count) {
    Frustum<Base> frustum(base, requested_height);
    Rect<T,1> h = frustum.clipped_height();
    Rect<T,3> bounds = frustum.bounds();
    P3 center = bounds.center();
    T radius = std::max<T>(bounds.dimensions().mag() / 2, 0.25);

    Sample sample {name, description, h, center, radius};
    add_wireframe(sample.wire, boundary, h);

    std::normal_distribution<T> normal;
    for (int i = 0; i < ray_count; ++i) {
        P3 origin = center + random_unit(rng) * radius * uniform(rng, 2.0, 4.0);
        P3 target = center + P3{normal(rng), normal(rng), normal(rng)} * radius * 0.75;
        P3 direction = (target - origin).unit();

        // Reserve a few lines for exact or nearly-degenerate branches; the
        // remainder are the randomized assault visible in the viewer.
        if (i == 0) {
            direction = -origin.unit(); // through the frustum tip
        } else if (i == 1) {
            origin[2] = center[2];
            direction = P3{random_unit(rng).template resized<2>().unit(), 0};
        } else if (i == 2) {
            origin = {2 * radius, 0, 0};
            direction = {-1, 0, 0}; // lies in h=0 and passes through the tip
        } else if (i == 3) {
            direction[2] *= 1e-10;
            direction = direction.unit();
        }
        Ray<T,3> ray {origin, direction};

        T s_center = (center - origin).dot(direction);
        Rect<T,1> draw_range {s_center - 3.5 * radius, s_center + 3.5 * radius};
        Rect<T,1> hit = frustum.intersect(ray) & draw_range;

        if (hit.is_empty()) {
            sample.outside.push_back({ray.at_multiple(draw_range.lo),
                                      ray.at_multiple(draw_range.hi)});
        } else {
            if (draw_range.lo < hit.lo) {
                sample.outside.push_back({ray.at_multiple(draw_range.lo),
                                          ray.at_multiple(hit.lo)});
            }
            sample.hits.push_back({ray.at_multiple(hit.lo), ray.at_multiple(hit.hi)});
            if (hit.hi < draw_range.hi) {
                sample.outside.push_back({ray.at_multiple(hit.hi),
                                          ray.at_multiple(draw_range.hi)});
            }
        }
    }
    return sample;
}

Sample random_rectangle(std::mt19937_64& rng, int rays) {
    P2 center {uniform(rng, -0.7, 0.7), uniform(rng, -0.7, 0.7)};
    P2 half {uniform(rng, 0.35, 1.5), uniform(rng, 0.35, 1.5)};
    Rect<T,2> base(center - half, center + half);
    std::vector<P2> boundary {
        {base.lo[0], base.lo[1]}, {base.hi[0], base.lo[1]},
        {base.hi[0], base.hi[1]}, {base.lo[0], base.hi[1]}
    };
    return make_sample("rectangle", "axis-aligned rectangular base", base, boundary,
                       random_height(rng), rng, rays);
}

Sample random_circle(std::mt19937_64& rng, int rays) {
    P2 center {uniform(rng, -0.7, 0.7), uniform(rng, -0.7, 0.7)};
    T radius = uniform(rng, 0.4, 1.5);
    Sphere<T,2> base(center, radius);
    std::vector<P2> boundary;
    constexpr int n = 48;
    for (int i = 0; i < n; ++i) {
        T a = 2 * std::numbers::pi_v<T> * i / n;
        boundary.push_back(center + radius * P2{std::cos(a), std::sin(a)});
    }
    return make_sample("circle", "circular base", base, boundary,
                       random_height(rng), rng, rays);
}

Sample random_triangle(std::mt19937_64& rng, int rays) {
    P2 center {uniform(rng, -0.5, 0.5), uniform(rng, -0.5, 0.5)};
    T angle = uniform(rng, 0, 2 * std::numbers::pi_v<T>);
    std::vector<P2> boundary;
    for (int i = 0; i < 3; ++i) {
        T a = angle + 2 * std::numbers::pi_v<T> * i / 3;
        T r = uniform(rng, 0.65, 1.5);
        boundary.push_back(center + r * P2{std::cos(a), std::sin(a)});
    }
    Simplex<T,2> base {boundary[0], boundary[1], boundary[2]};
    return make_sample("triangle", "triangular simplex base", base, boundary,
                       random_height(rng), rng, rays);
}

void write_point(std::ostream& out, const P3& p) {
    out << '[' << p[0] << ',' << p[1] << ',' << p[2] << ']';
}

void write_segments(std::ostream& out, const std::vector<Segment>& segments) {
    out << '[';
    for (std::size_t i = 0; i < segments.size(); ++i) {
        if (i) out << ',';
        out << '[';
        write_point(out, segments[i].first);
        out << ',';
        write_point(out, segments[i].second);
        out << ']';
    }
    out << ']';
}

void write_samples(const Options& options, const std::vector<Sample>& samples) {
    std::ofstream out(options.output);
    if (!out) throw std::runtime_error("could not open output file: " + options.output);
    out << std::setprecision(17);
    out << "{\"seed\":" << options.seed << ",\"rays_per_sample\":" << options.rays
        << ",\"samples\":[";
    for (std::size_t i = 0; i < samples.size(); ++i) {
        if (i) out << ',';
        const Sample& s = samples[i];
        out << "{\"shape\":\"" << s.shape << "\",\"description\":\""
            << s.description << "\",\"height\":[" << s.height.lo << ',' << s.height.hi
            << "],\"center\":";
        write_point(out, s.center);
        out << ",\"radius\":" << s.radius << ",\"wire\":";
        write_segments(out, s.wire);
        out << ",\"outside\":";
        write_segments(out, s.outside);
        out << ",\"hits\":";
        write_segments(out, s.hits);
        out << '}';
    }
    out << "]}\n";
}

std::vector<std::string> split_shapes(std::string_view text) {
    std::vector<std::string> result;
    std::size_t start = 0;
    while (start <= text.size()) {
        std::size_t comma = text.find(',', start);
        std::string value(text.substr(start, comma == std::string_view::npos
            ? text.size() - start : comma - start));
        if (!value.empty()) result.push_back(value);
        if (comma == std::string_view::npos) break;
        start = comma + 1;
    }
    return result;
}

Options parse_options(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (++i == argc) throw std::runtime_error("missing value for " + arg);
            return argv[i];
        };
        if (arg == "--output") o.output = value();
        else if (arg == "--samples") o.samples = std::stoi(value());
        else if (arg == "--rays") o.rays = std::stoi(value());
        else if (arg == "--seed") o.seed = std::stoull(value(), nullptr, 0);
        else if (arg == "--shapes") o.shapes = split_shapes(value());
        else if (arg == "--help") {
            std::cout << "Usage: frustum_visualizer_data [options]\n"
                      << "  --output PATH\n  --samples N\n  --rays N\n"
                      << "  --seed N\n  --shapes rectangle,circle,triangle\n";
            std::exit(0);
        } else throw std::runtime_error("unknown option: " + arg);
    }
    if (o.samples < 1 || o.rays < 1 || o.shapes.empty()) {
        throw std::runtime_error("samples, rays, and shape selection must be nonempty");
    }
    for (const std::string& shape : o.shapes) {
        if (shape != "rectangle" && shape != "circle" && shape != "triangle") {
            throw std::runtime_error("unknown shape: " + shape);
        }
    }
    return o;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options options = parse_options(argc, argv);
        std::mt19937_64 rng(options.seed);
        std::vector<Sample> samples;
        samples.reserve(options.samples);
        for (int i = 0; i < options.samples; ++i) {
            const std::string& shape = options.shapes[i % options.shapes.size()];
            if (shape == "rectangle") samples.push_back(random_rectangle(rng, options.rays));
            else if (shape == "circle") samples.push_back(random_circle(rng, options.rays));
            else samples.push_back(random_triangle(rng, options.rays));
        }
        std::shuffle(samples.begin(), samples.end(), rng);
        write_samples(options, samples);
        std::cout << "Wrote " << samples.size() << " samples to " << options.output
                  << " (seed " << options.seed << ")\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
