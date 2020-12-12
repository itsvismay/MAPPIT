#pragma once

#include <Eigen/Sparse>

#include <cmath>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <set>
#include <vector>

#include "bvh.h"

double smoothExpDist(double d, double alpha, double dmax);

std::pair<int, double> findClosestPoint(
    const Eigen::VectorXd& p, int idx, const BVH* bvh,
    double maxRadius = std::numeric_limits<double>::infinity());

struct SmoothDistComponent {
  double dist;
  Eigen::VectorXd grad;

  SmoothDistComponent(double dist, const Eigen::VectorXd& grad)
      : dist(dist), grad(grad) {}

  SmoothDistComponent()
      : dist(0), grad(Eigen::VectorXd::Zero(3)) {}

  SmoothDistComponent operator+(const SmoothDistComponent& b) {
    SmoothDistComponent c = *this;
    c.dist += b.dist;
    c.grad += b.grad;
    return c;
  }

  SmoothDistComponent& operator+=(const SmoothDistComponent& b) {
    dist += b.dist;
    grad += b.grad;
    return *this;
  }
};

struct SmoothDistResult {
  double alpha;
  double dmax;
  double true_dist;
  double smooth_dist;
  Eigen::VectorXd grad; // wrt query point
  int retained_pairs;
};

double goodAlpha(double true_min, double c_upper, double d_max);

SmoothDistComponent smoothPointExpDist(const Eigen::VectorXd& p, int idx,
                                       const BVH* bvh, double alpha,
                                       double dmax, std::vector<int>& num_used);

SmoothDistResult smoothMinDist(const Eigen::MatrixXd& V, const BVH* bvh, double c_upper, double box_frac, const Eigen::VectorXd& p);
