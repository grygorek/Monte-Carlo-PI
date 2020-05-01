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
// This example implements two methods of estimating PI.
//  - estimating an area of the circle
//  - averaging value of a circle function


// Area Method
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
auto pi_area_method()
{
  std::random_device rnd_dev;
  std::mt19937 gen(rnd_dev());
  std::uniform_real_distribution<double> rnd(0, std::nextafter(1.f, 2.f));

  auto area_method = [&gen, &rnd](int nloops) {
    int inside{};
    while (nloops--)
    {
      auto point = std::complex<double>{rnd(gen), rnd(gen)};
      if (std::abs(point) <= 1.0)
        inside++;
    }
    return inside;
  };

  constexpr int n = 10;
  constexpr int all_points = 10'000;

  // just for fun spin few threads
  std::vector<std::future<int>> tasks(n);
  for (auto& t : tasks)
    t = std::async(area_method, all_points / n);

  int inside = std::accumulate(
      std::begin(tasks), std::end(tasks), 0, [](auto acc, auto& i) {
        return acc + i.get();
      });

  return 4. * inside / all_points;
}

// Averaging value of a circle function method
//
// The circle function is r^2 = x^2 + y^2. If radius r is 1 then:
//   y = sqrt(1 - x^2)
//
// The integral of that function is Pi/4.
//  Favg = 1/(b-a) * SUM[i=a..b](sqrt(1 - Xi^2))
// In our case a=0 and b=1 so:
//  Favg = SUM[i=0..1](sqrt(1 = Xi^2))
// In other words,
// generate lots of random x, sum up all y and divide by number of samples.
auto pi_avg_method()
{
  std::random_device rnd_dev;
  std::mt19937 gen(rnd_dev());
  std::uniform_real_distribution<double> rnd(0, std::nextafter(1.f, 2.f));

  auto avg_method = [&gen, &rnd](int nloops) {
    double y{};
    while (nloops--)
    {
      auto x = rnd(gen);
      y += sqrt(1 - x * x);
    }
    return y;
  };

  constexpr int n = 10;
  constexpr int all_points = 10'000;

  // Just for fun, run concurently.
  std::vector<std::future<double>> tasks(n);
  for (auto& t : tasks)
    t = std::async(avg_method, all_points / n);

  auto y = std::accumulate(
      std::begin(tasks), std::end(tasks), 0., [](auto acc, auto& i) {
        return acc + i.get();
      });

  return 4. * y / all_points;
}

template <class T, class Type>
void histogram(const T& dta, Type min, Type max, int bins_cnt, int scale = 1)
{
  // For each value find the correct bin and increment its count.

  std::vector<int> bins(bins_cnt);
  const auto step = (max - min) / bins_cnt;
  for (auto v : dta)
  {
    for (auto bin = 0; bin < bins_cnt; bin++)
    {
      auto mn = min + step * bin;
      auto mx = min + step * (bin + 1LL);
      if (v >= mn && v < mx)
      {
        ++bins[bin];
        break;
      }
    }
  }

  // print min/max of the histogram (left and right edge)
  std::cout << "min: " << min << ", max: " << max << '\n';
  
  // for each bin print a bar as long as the bin's value
  // (use scale argument if line in the console is broken)
  for (int bin = 0; bin < bins_cnt; bin++)
    std::cout << bin << ": " << std::string(bins[bin] / scale, '*') << '\n';
}

template <class F>
void do_statistics_pi(F&& f)
{
  // Generate PI N times and print mean value, standard error
  // and histogram.
    
  std::vector<double> pi_est(500);

  double sum2{};  // sum of squares (needed for std error)
  double sum{};
  for (auto& pi : pi_est)
  {
    pi = f();
    sum += pi;
    sum2 += pi * pi;
  }

  // To speed up the variation calculation (= std error) we can use a property
  // that variation is a difference between a avg of squered samples and a
  // squered sum of means var = SUM( sum^2 ) - SUM( mean )^2

  sum2 = sum2 / pi_est.size();  // avg of squared samples
  double u = sum / pi_est.size();  // mean
  double var = sum2 - u * u;  // avg of sqr smpls - mean^2
  double std_err = sqrt(var);

  // This is a slow version of variation (loop around samples again is needed).
  // Doing this for comparison as a proof there is no difference.
  double v{};
  for (auto pi : pi_est)
  {
    auto t = pi - u;
    v += t * t;
  }
  v /= pi_est.size();
  auto std_err2 = sqrt(v);

  std::cout << "Estimated PI: " << u << ", std err(fast): " << std_err
            << ", std err2(slow): " << std_err2 << "\n";

  // print a histogram
  const auto [min, max] = std::minmax_element(pi_est.begin(), pi_est.end());
  histogram(pi_est, *min, *max, 10, 5);
}

int main()
{
  std::cout << "Area Method:\n";
  do_statistics_pi(pi_area_method);

  std::cout << "\nAverage Method:\n";
  do_statistics_pi(pi_avg_method);
}
