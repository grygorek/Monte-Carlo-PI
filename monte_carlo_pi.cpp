#include <complex>
#include <future>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

// Estimating PI number using Monte Carlo method
//
// Bit of maths.
//
// 1. What is PI?
//
// It is a ratio between a circumference (C) of a circle and a diameter (d) of
// that circle. Pi = C / d
//
// Also it can be found from the circle area equation and a given radious (r):
// Pi = Area / (r*r) Geometrically, r*r is an area of a square with a side
// length of 'r'. So, PI is a ratio between a circle area to the square area.
//
// If that square is drawn over a circle, one corner of the square being in the
// centre of the circle, then the other two corners will overlap with the edges
// of the circle. It is easy to notice that the area of the circle overlapping
// with that square is exactly 1/4 of the circle. Or... you need 4 such squares
// to cover the entire circle. Either way, the ratio is 4.
//
// If, CircArea4 is circle area within a single square, then:
//  PI = 4 * CircArea4 / square_area
//
// That means that to calculate the PI the ratio of the areas is needed.
//
// 2. What is Monte Carlo metod?
//
// Very simply, it is a numerical metod to estimate results using random
// sampling.
//
// 3. Using Monte Carlo method to estimate the ratio between circle and square
// areas
//
// Say, the radious of the circle is 1.0. Consequently, the square is 1x1. Let's
// generate random points on the surface of that square. Some of them will also
// overlap with the surface of the quarter of the circle and the others will
// fall outside of the circle. Points within the circle have the distance from
// the square corner (circle centre) less than 1.0.
//
// The ratio of the number of points falling on the circle (1/4 of the total
// circle area) to the number of all points (square area) is the ratio we are
// looking for.
//
// ratio = circle_points / all_points
// is the same as
// ratio = CircArea4 / square_area
//
// Finally, estimated PI = 4 * ratio

int main() {
  std::random_device rnd_dev;
  std::mt19937 gen(rnd_dev());
  std::uniform_real_distribution<double> rnd;

  auto calculate = [&gen, &rnd](int nloops) {
    int inside{};
    while (nloops--) {
      auto point = std::complex<double>{rnd(gen), rnd(gen)};
      if (std::abs(point) <= 1.0) inside++;
    }
    return inside;
  };

  constexpr int n = 10;
  constexpr int all_points = 5'000'000;

  std::vector<std::future<int>> tasks(n);
  for (auto& t : tasks) t = std::async(calculate, all_points / n);

  int inside = std::accumulate(std::begin(tasks), std::end(tasks), 0,
                               [](auto acc, auto& i) { return acc + i.get(); });

  std::cout << "Estimated PI: " << 4.0 * inside / all_points << "\n";
}
