


//  Purpose: Detecting turning points in high-frequency animal movement data
// 
//  Author: Abdulmajeed Alharbi
// 
//  Reference: ....
// 
//  Date created: January 11 2025



#include <Rcpp.h>
#include <boost/numeric/interval.hpp>
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_set.hpp>
#include <boost/math/constants/constants.hpp>

using namespace Rcpp;
using namespace boost::numeric;
using namespace boost::icl;




// Create a closed interval
interval_set<double> CreateInterval(double s, double e) {
  interval_set<double> result;
  result += boost::icl::interval<double>::closed(s, e);
  return result;
}

// Calculate the complement of an interval set in the circular range [-π, π]
interval_set<double> complementInterval(interval_set<double> set) {
  interval_set<double> range = CreateInterval( -M_PI, M_PI );
  interval_set<double> result = range - set;
  return result;
}


// Class to handle and store cosine wave segments
// Each cosine wave is characterised by 4 components 
// R cos(mu - alpha) +  z      in the range I

class CosineSegment {

    public:
      CosineSegment(interval_set<double> I , int tau, double a, double b, double c): I(I), tau(tau), a(a), b(b), c(c) {
        R = 1;                        // Initialize resultant length
        maxima = R + c;               // Initialize maximum based on R
      }
      
      // Update the segment with new sine and cosine values
      void update(double sin_h, double cos_h) {
        a += sin_h;
        b += cos_h;
        R = std::sqrt(a * a + b * b);
        maxima = R + c;
      }
      
      // Find the solution intervals satisfying a penalized value
      interval_set<double> root_finder(double z) {
        double discriminant = (z - c) / R; 
        
        // If the value is out of bounds, return empty interval (No solution)
        if (discriminant > 1) {
          return interval_set<double>(); // Return an empty interval
        }
        
        // If the solution spans the entire space.
        if (discriminant < -1) {
          return I;
        }
        
        // Otherwise, compute the solution interval
        double alpha = std::atan2(a, b);    // Mean angle
        double lhs   = std::acos(discriminant);      // Angular displacement
        
        // Calculate the lower and upper bounds for the interval
        double lower = std::fmod(-lhs + alpha + M_PI, 2 * M_PI) - M_PI;
        double upper = std::fmod(lhs  + alpha + M_PI, 2 * M_PI) - M_PI;
        
        interval_set<double> result;
        if (lower > upper) {
          interval_set<double> A = CreateInterval(upper, lower);
          result = complementInterval(A);
        } else {
          result = CreateInterval(lower, upper);
        }
        
        I =  result&I;
        return I;
      }
      
      interval_set<double> I;  // Interval for the segment
      int tau;                 // Time index
      double a;                // Sine component
      double b;                // Cosine component
      double c;                // Constant value
      double R;                // Resultant length (magnitude of vector)
      double maxima;           // Maximum value (R + c)
};








// [[Rcpp::export]]
std::vector<int> AFPOP(NumericVector sines, NumericVector cosines, double penalty) {
  // Parameters:
  // sines - Sine values of headings
  // cosines - Cosine values of headings
  // penalty - Penalty constant for segmentation
  
  
  int n = sines.size();
  interval_set<double> Range = CreateInterval(-M_PI, M_PI);
  CosineSegment F1(Range, 0, sines[0], cosines[0], penalty);
  
  // To store the details for each segment
  std::vector<CosineSegment> Fs; 
  Fs.emplace_back(F1);
  
  // A vector to store the most recent turning-points
  std::vector<int> last_tps(n);
    
  for(int t = 1; t < n; ++t) {
    double maximum = -std::numeric_limits<double>::infinity();
    int tau = 0;
    int nn = Fs.size();
    
    // To calculate F(t) and the most probable turing-points
    for(int j = 0; j < nn; ++j) {
      if(Fs[j].maxima > maximum) {
        maximum = Fs[j].maxima;
        tau = Fs[j].tau;
      }
    }
    
    
    
    interval_set<double> Union;
    std::vector<CosineSegment> Fs_new;
    last_tps[t] = tau;
    maximum += penalty;
    
    for(int j = 0; j < nn; ++j) {
      interval_set<double> root = Fs[j].root_finder(maximum);
      if(!root.empty()) {
        Union = (Union | root);
        Fs[j].update(sines[t], cosines[t]);
        Fs_new.emplace_back(Fs[j]);
      }// else purn it. 
    }
    
    interval_set<double> new_range = complementInterval(Union);
    CosineSegment Ft(new_range, t, sines[t], cosines[t], maximum);
    Fs_new.emplace_back(Ft);
    Fs = Fs_new;
    
  }
  
  
  std::vector<int> tps = {n}; // Initialize tps with n
  int last = n-1;
  while(last>0){
    last = last_tps[last];
    tps.push_back(last);
  }
  
  return tps;
}










// [[Rcpp::export]]
std::vector<int> DFPOP(const NumericMatrix& M,
                       const std::vector<double>& D,
                       const double rho,
                       const double penalty) {
  
  // The implementation of AR1 FPOP algorithm via discretisation
  
  int d = M.nrow();    // Number of rows in M
  int n = M.ncol();    // Number of columns in M
  
  std::vector<double> likelihood_curve(d, 0.0);  
  std::vector<int> range(d, 0);
  std::vector<int> last_tps(n, 0);
  double maximum = penalty;
  double q = rho/(1-rho);
  double qq = 1 + q * q;
  double q2 = 2 * q;
  
  for (int t = 1; t < n; ++t) {
    double maxx = -std::numeric_limits<double>::infinity();
    int argmaxx = -1;
    
    for (int j = 0; j < d; ++j) {
      double numerator = M(j,t) + q * D[t - 1];
      double denominator =  sqrt(qq + q2 * M(j, t-1));
      double cosine_curve = numerator / denominator;
      
      if (likelihood_curve[j] <= maximum) {
        likelihood_curve[j] = maximum;
        range[j] = t;
      }
      
      likelihood_curve[j] += cosine_curve;
      
      if (likelihood_curve[j] > maxx) {
        maxx = likelihood_curve[j];
        argmaxx = j;
      }
    }
    
    maximum = penalty + maxx;
    last_tps[t] = range[argmaxx];
  }
  
  
  std::vector<int> tps = {n};
  int last = n-1;
  while(last>1){
    last = last_tps[last-1];
    tps.push_back(last);
  }
  
  return tps;
}

