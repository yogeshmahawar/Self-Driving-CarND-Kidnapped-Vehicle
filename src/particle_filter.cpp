/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

#define NUMBER_OF_PARTICLES 20 // Can be decreased (even 12 particles can pass the test)
#define EPS 0.001  // Just a small number

void ParticleFilter::init(double x, double y, double theta, double std[]) {

	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	//generate random number for initializing particle x,y and theta in range of gps std values
	static default_random_engine gen;
	//creating normal distribution of initial gps values
	// Add random Gaussian noise to each particle.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	//setting up number of particles
	num_particles = NUMBER_OF_PARTICLES;
	for(unsigned short  int i=0 ; i<num_particles ; i++){
			Particle  particle = Particle();
			particle.x = dist_x(gen);
			particle.y = dist_y(gen);
			particle.theta = dist_theta(gen);
			//to be checked
			particle.weight = 0.0;
			particles.push_back(particle);
	}
  is_initialized =true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

	const double vdt = velocity * delta_t;
	const double yawdt = yaw_rate * delta_t;
	const double vel_yaw = velocity/yaw_rate;
	normal_distribution<double> dist_x(0.0, std_pos[0]);
  normal_distribution<double> dist_y(0.0, std_pos[1]);
  normal_distribution<double> dist_theta(0.0, std_pos[2]);
	static default_random_engine gen;
	for (unsigned short int i = 0; i < num_particles; i++){
		    if (fabs(yaw_rate) < EPS){
            particles[i].x += vdt * cos(particles[i].theta);
            particles[i].y += vdt * sin(particles[i].theta);
            // particles[i].theta unchanged if yaw_rate is too small
						//fomula is bit diffeent as below
						// xf = x0 + vdt + cos(theta)
						// yf = y0 + vdt + cos(theta)
        }
        else{
            const double theta_new = particles[i].theta + yawdt;
            particles[i].x += vel_yaw * (sin(theta_new) - sin(particles[i].theta));
            particles[i].y += vel_yaw * (-cos(theta_new) + cos(particles[i].theta));
            particles[i].theta = theta_new;
						//when yaw_rate is ok, then
						//xf = x0 + (vel/yaw_rate)( sin(theta_new) - sin(old_theta))
						//yf = y0 + (vel/yaw_rate)( -cos(theta_new) + cos(old_theta))
						//theta = new_theta
        }
        // Add random Gaussian noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	//calculating sigma2 for both x and y positions
	const double sigma_xx = std_landmark[0] * std_landmark[0];
	const double sigma_yy = std_landmark[1] * std_landmark[1];
	//calculating normalization factor for multivariate gaussian
	const double normal_factor = 2 * M_PI * std_landmark[0] * std_landmark[1];
	double dx = 0.0;
	double dy = 0.0;
	double sum_w = 0.0;
	weights.resize(num_particles);
	//iterate through all particles
	for(unsigned short int i=0 ; i< num_particles ; i++){
		unsigned short int counter = 0;
		double weight = 0.0;
		const double sin_theta = sin(particles[i].theta);
		const double cos_theta = cos(particles[i].theta);
		//iterating through all observations
		for (unsigned short int j = 0; j < observations.size(); j++){
			// Observation measurement transformations with formula equations below
			LandmarkObs obs;
			obs.id = observations[j].id;
			// xm = xp + cos(theta)*xc - sin(theta)*yc
			obs.x = particles[i].x + (observations[j].x * cos_theta) - (observations[j].y * sin_theta);
			// ym = yp + sin(theta)*xc + cos(theta)*yc
			obs.y = particles[i].y + (observations[j].x * sin_theta) + (observations[j].y * cos_theta);
			// Unefficient way for observation asossiation to landmarks. It can be improved.
			Map::single_landmark_s nearest_lm;
			double nearest_dist = 10000000.0; // A big number
			//finding out nearest land mark
      for (unsigned short int k = 0; k < map_landmarks.landmark_list.size(); k++) {
				Map::single_landmark_s cond_lm = map_landmarks.landmark_list[k];
				//calculating distance between landmark and particle observation
				// to idetify nearest landmark to the observation
        double distance = dist(cond_lm.x_f, cond_lm.y_f, obs.x, obs.y);  // Calculate the Euclidean distance between two 2D points
        if (distance < nearest_dist) {
					nearest_dist = distance;
					// to idetify nearest landmark to the observation
          nearest_lm = cond_lm;
          // if (distance < sensor_range){
					// 	in_range = true;
					// }
        }
		  }
			double distance_from_particle = dist(nearest_lm.x_f, nearest_lm.y_f, particles[i].x, particles[i].y);
			if (distance_from_particle <= sensor_range){
				//nearest and also in range
				dx = obs.x-nearest_lm.x_f;
				dy = obs.y-nearest_lm.y_f;
				weight += ((dx * dx )/ sigma_xx) + ((dy * dy) / sigma_yy);
				counter++;
			}

		}
		//particles[i].weight = pow(normal_factor,counter) * exp(-0.5*weight); // calculate exp() after main computation in order to optimize the code
		//removing normal factor because it is getting devided in normalization step
		particles[i].weight = exp(-0.5*weight);
		//cout << "weight of particles : "  << particles[i].weight << endl;
		weights[i] = particles[i].weight;
		//cout << "weight of particle :" << weights[i];
    sum_w += weights[i];
  }
	for (unsigned short int i = 0; i < num_particles; i++){
		particles[i].weight /= sum_w;
		weights[i] = particles[i].weight;
		//cout << "weight of particles : "  << particles[i].weight << endl;
	}

}

void ParticleFilter::resample() {
	vector<Particle> new_particles;
  new_particles.resize(num_particles);
  double max_w = *std::max_element( weights.begin(), weights.end());
	//cout << "max weight in resampling " << max_w;
	unsigned short int index = std::rand() % (num_particles);
	//cout << "random index : " << index;
	double beta = 0.0 ;
	for (unsigned short int i = 0; i<num_particles; i++){
		beta = beta + 2.0 * max_w;
    while (beta > weights[index]){
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
    new_particles[i] = particles[index];
	}
	//particles = new_particles;
	// for (int i = 0; i < num_particles; i++){
	// 	particles[i] = new_particles[i];
	// 	//cout << "weight of particles after resampling : "  << particles[i].weight << endl;
	// }
	particles.assign(new_particles.begin(), new_particles.end());
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
