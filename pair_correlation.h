#ifndef GUARD_PAIR_CORRELATION_H
#define GUARD_PAIR_CORRELATION_H


#include <vector>
#include <math.h>
#include <fstream>
#include <string>

#include "vec3.h"


class PairCorrelation {

public:
	// constructor
	PairCorrelation(unsigned int number_of_bins, double bins_size,
					double bulk_density, double system_size_x,
					double system_size_y, double system_size_z);

	// sample
	void sample( const std::vector<Vec3>& particles );

    std::vector<double> GetPairCorrelation() const;
    std::vector<double> GetBins() const {return bins;}

	// get current number of samples taken
	long unsigned int GetNumberOfSamples() { return number_of_samples_;}
	long unsigned int GetNumberOfSamplesIgnored() { return number_of_samples_ignored_;}

	// write to out stream
	void write(std::ostream & out);
	// write to file called outname
	void write(std::string outname);
private:

	unsigned int number_of_bins_;
	double bin_size_;
  double bulk_density_;

	double system_size_x_;
	double system_size_y_;
	double system_size_z_;


  double max_distance_; // if distance is larger than this
                        // does not contribute to the pair 
                        // correlation

	long unsigned int number_of_samples_;	// number of samples

  // number of samples with distance > max_distance
  long unsigned int number_of_samples_ignored_;

	std::vector<double> bins;	// mid values of the bins 
	std::vector<double> gr;		// pair correlation function
	
};


PairCorrelation::PairCorrelation(
		unsigned int number_of_bins,
		double bin_size,
    double bulk_density,
		double system_size_x,
		double system_size_y,
		double system_size_z)
  : number_of_bins_(number_of_bins),
    bin_size_(bin_size),
    bulk_density_(bulk_density),
    system_size_x_(system_size_x),
    system_size_y_(system_size_y),
    system_size_z_(system_size_z),
    max_distance_(bin_size * number_of_bins),
    number_of_samples_(0),
    number_of_samples_ignored_(0),
    bins(number_of_bins),
    gr(number_of_bins, 0.0)
{

	// set mid points of the bins
	for(unsigned int i = 0; i < number_of_bins_; ++i)
		bins[i] = (0.5 + i) * bin_size_;
}			

void PairCorrelation::sample( const std::vector<Vec3>& particles )
{
	for(unsigned int i = 0; i < particles.size(); ++i) {
		for(unsigned int j = i + 1; j < particles.size(); ++j) {
      Vec3 r = particles[i] - particles[j]; 
      // periodic boundary conditions
      if (system_size_x_ > 0) {
        r.x -= system_size_x_ * round(r.x / system_size_x_);
      }
      if (system_size_y_ > 0) {
        r.y -= system_size_y_ * round(r.y / system_size_y_);
      }
      if (system_size_z_ > 0) {
        r.z -= system_size_z_ * round(r.z / system_size_z_);
      }
      double d = r.Length();
      if(d < max_distance_) {
        gr[ (int)(d / bin_size_)] += 2.0;
        number_of_samples_ += 1;
      } else {
        number_of_samples_ignored_ += 1;
      }
    }
  }
}

std::vector<double> PairCorrelation::GetPairCorrelation() const
{
    std::vector<double> gr_current = gr;

	const double pi = 4 * atan(1.0);
	// number of ideal gas particles in bin.
	//double ideal_number_of_particles;
  double bulk_density = system_size_x_ * system_size_y_ * system_size_z_; 
  for (unsigned int i = 0; i < number_of_bins_; ++i) {
    double volume = pow( (1 + i) * bin_size_, 3.0);
    volume -= pow( i * bin_size_, 3.0);
    volume *= 4.0 * pi / 3.0;
    gr_current[i] = gr[i];
    gr_current[i] /= number_of_samples_ * volume * 2 / bulk_density;
  }
  
	//for(unsigned int i = 0; i < number_of_bins_; ++i) {
  //      double bin_volume = pow((1 + i) * bin_size_, 3.0) - pow(i * bin_size_, 3.0);
  //      bin_volume *= 4.0 * pi / 3.0;

  //      ideal_number_of_particles = bin_volume * bulk_density_;

	//	gr_current[i] = gr[i];
	//	gr_current[i] /= ideal_number_of_particles * number_of_samples_; 
	//}	
  return gr_current;
}


void PairCorrelation::write(std::ostream & out)
{
  std::vector<double> gr_current = GetPairCorrelation();

	for(unsigned int i=0;i < number_of_bins_;++i) {
		out << bins[i] << ' ';
		out << gr_current[i] << '\n';
	}

}


void PairCorrelation::write(std::string outname)
{
  std::vector<double> gr_current = GetPairCorrelation();

	std::ofstream out;
	out.open(outname);
	for(unsigned int i = 0; i < number_of_bins_; ++i) {
		out << bins[i] << ' ' << gr_current[i];
		if( i < number_of_bins_ - 1) out <<  '\n';
	}
	out.close();
}
#endif
