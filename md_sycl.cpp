/* A serial example molecular dynamics program
*
* Tested using hipSYCL
*
* compile using
* syclcc -lm md_sycl.cpp -o md_sycl
*
* run using
* ./md_sycl
*/
#include<CL/sycl.hpp>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

void forceupdate(float *position, float *force, 
                int numparticles, float sigma, float epsilon)
{

    // start device queue
    cl::sycl::default_selector deviceSelector;
    cl::sycl::queue queue(deviceSelector);
    // Create device memory
    cl::sycl::buffer<float, 1> d_position(position , numparticles*3);
    cl::sycl::buffer<float, 1> d_force(force , numparticles*3);
    queue.submit([&d_position, &d_force, sigma, epsilon, numparticles](
    cl::sycl::handler& cgh) {
      auto force_data = d_force.get_access
      <cl::sycl::access::mode::write>(cgh);
     cgh.parallel_for(cl::sycl::range<1>(3*numparticles),
        [force_data](cl::sycl::id<1> idx){
        force_data[idx] = 0;
        });
     });
    queue.wait();
    queue.submit([&d_position, &d_force, sigma, epsilon, numparticles](
    cl::sycl::handler& cgh) {
      auto force_data = d_force.get_access
      <cl::sycl::access::mode::read_write>(cgh);
     auto position_data = d_position.get_access
      <cl::sycl::access::mode::read>(cgh);
      cgh.parallel_for<class force>(
        cl::sycl::nd_range<1>{numparticles,numparticles},
        [position_data, force_data, sigma, epsilon, numparticles](cl::sycl::nd_item<1> item){
          int global_id = item.get_global_linear_id();
          int ii = global_id % numparticles;
          int jj = (global_id - ii) / numparticles;
          float rsq = (position_data[3*ii  ]-position_data[3*jj  ])*
                      (position_data[3*ii  ]-position_data[3*jj  ])
                     +(position_data[3*ii+1]-position_data[3*jj+1])*
                      (position_data[3*ii+1]-position_data[3*jj+1])
                     +(position_data[3*ii+2]-position_data[3*jj+2])*
                      (position_data[3*ii+2]-position_data[3*jj+2]);
       float mindist = 0.000000001;
       rsq = ( rsq > mindist ) ? rsq : mindist;   
      // x,y and z components
        force_data[3*ii  ] +=
             24 * epsilon * ( (position_data[3*ii  ]-position_data[3*jj  ]) / (sigma*sigma) )  *
            ( 2 * pow(( sigma*sigma / rsq ),7) -  pow(( sigma*sigma / rsq ),4) );
        force_data[3*ii+1] +=
             24 * epsilon * ( (position_data[3*ii+1]-position_data[3*jj+1]) / (sigma*sigma) )  *
            ( 2 * pow(( sigma*sigma / rsq ),7) -  pow(( sigma*sigma / rsq ),4) );
        force_data[3*ii+2] +=
             24 * epsilon * ( (position_data[3*ii+2]-position_data[3*jj+2]) / (sigma*sigma) )  *
            ( 2 * pow(( sigma*sigma / rsq ),7) -  pow(( sigma*sigma / rsq ),4) );
     });
    });
    queue.wait();  
 
  return;
}

void positionupdate(float *position, float *velocity, 
                    float dt, int numparticles)
{
    // start device queue
    cl::sycl::default_selector deviceSelector;
    cl::sycl::queue queue(deviceSelector);
    // Create device memory
    cl::sycl::buffer<float, 1> d_position(position , numparticles*3);
    cl::sycl::buffer<float, 1> d_velocity(velocity , numparticles*3);
    queue.submit([&d_position, &d_velocity, dt, numparticles](
    cl::sycl::handler& cgh) {
      auto position_data = d_position.get_access
      <cl::sycl::access::mode::read_write>(cgh);
     auto velocity_data = d_velocity.get_access
      <cl::sycl::access::mode::read>(cgh);
      cgh.parallel_for(cl::sycl::range<1>(3*numparticles),
        [position_data, velocity_data, dt](cl::sycl::id<1> idx){
        position_data[idx] += dt * velocity_data[idx];
        });
      });
    queue.wait();
   return;
}

void velocityupdate(float *velocity, float *force,
                    float dt, float mass, int numparticles)
{
    // start device queue
    cl::sycl::default_selector deviceSelector;
    cl::sycl::queue queue(deviceSelector);
    // Create device memory
    cl::sycl::buffer<float, 1> d_force(force , numparticles*3);
    cl::sycl::buffer<float, 1> d_velocity(velocity , numparticles*3);
    queue.submit([&d_force, &d_velocity, dt,mass, numparticles](
    cl::sycl::handler& cgh) {
      auto force_data = d_force.get_access
      <cl::sycl::access::mode::read>(cgh);
     auto velocity_data = d_velocity.get_access
      <cl::sycl::access::mode::read_write>(cgh);
      cgh.parallel_for(cl::sycl::range<1>(3*numparticles),
        [force_data, velocity_data, dt, mass](cl::sycl::id<1> idx){
        velocity_data[idx] += dt * force_data[idx] / mass;
        });
      });
    queue.wait();
    return;
}

float kineticenergy(float *velocity, float mass, int numparticles)
{
  float ek=0.0;
  for(int i=0;i<3*numparticles;i++)
  {
    ek+=velocity[i]*velocity[i];
  }
  ek=0.5*mass*ek;
  return ek;
}

float potentialenergy(float *position, float epsilon, float sigma, int numparticles)
{
        float ep=0.0;
        for(int i=0;i<numparticles;i++)
        {
          for(int j=0;j<i;j++)
          {
            float rsq = (position[3*i  ]-position[3*j  ])*
                        (position[3*i  ]-position[3*j  ])
                      + (position[3*i+1]-position[3*j+1])*
                        (position[3*i+1]-position[3*j+1])
                      + (position[3*i+2]-position[3*j+2])*
                        (position[3*i+2]-position[3*j+2]); 
          
             ep+=pow(sigma*sigma/rsq,6) - pow(sigma*sigma/rsq,3);
          }
          for(int j=i+1;j<numparticles;j++)
          {
            float rsq = (position[3*i  ]-position[3*j  ])*
                        (position[3*i  ]-position[3*j  ])
                      + (position[3*i+1]-position[3*j+1])*
                        (position[3*i+1]-position[3*j+1])
                      + (position[3*i+2]-position[3*j+2])*
                        (position[3*i+2]-position[3*j+2]);

             ep+=pow(sigma*sigma/rsq,6) - pow(sigma*sigma/rsq,3);
          }
        }
        ep=2*epsilon*ep;
        return ep;
}

int main()
{
  // declare program parameters
  const int numparticles=50;
  const float epsilon=1;
  const float sigma=0.1;
  const int timesteps=10000;
  const float dt=0.001;
  const float mass=1.0;
  const int randinit=431;
  const int randmod=32767;
  const int enerdispsteps=100;
  // Declare arrays
  float position[3*numparticles];
  float velocity[3*numparticles];
  float force[3*numparticles];
  float ek;
  float ep;
  float et;
  time_t start;
  time_t stop;
  
  // Initialize random number generator
  srand(randinit);
  // Set initialize positions randomly
  for(int i=0; i < 3*numparticles; i++)
  {
    float temp = 0.3 + 0.01*i*i +
            0.0001*((float) (rand() % randmod)) / ((float) (randmod));
    position[i]=temp;
  } 

  // Zero initial velocity
  for(int i = 0; i < 3*numparticles; i++) 
  {
    velocity[i]=0;
  }
 
   // Calculate potential energy
   ep = potentialenergy(position, epsilon, sigma, numparticles);
   // Calculate kinetic energy
   ek = kineticenergy(velocity, mass, numparticles);
   // Get total energy
   et = ek+ep;
   printf("Time %f of %f \n Kinetic energy %f \n Potential Energy %f \n Total energy %f \n",
                    0*dt,dt*timesteps,ek,ep,et);
   // Update forces 
   forceupdate(position, force, numparticles, sigma, epsilon);
  start=time(NULL);
   // Time stepping loop
   for(int n = 0; n < timesteps; n++)
   {
     // Calculate new velocities 
     velocityupdate(velocity, force, 0.5*dt, mass, numparticles);
     // Calculate new positions
     positionupdate(position, velocity, dt, numparticles);
     // Update forces 
     forceupdate(position, force, numparticles, sigma, epsilon);
     // Calculate new velocities
     velocityupdate(velocity, force, 0.5*dt, mass, numparticles);

     if(n%enerdispsteps == 0)
     {
      // Calculate potential energy
      ep = potentialenergy(position,epsilon, sigma, numparticles);
      // Calculate kinetic energy
      ek = kineticenergy(velocity, mass, numparticles);
      // Get total energy
      et = ep+ek;
      printf("Time %f of %f \n Kinetic energy %f \n Potential Energy %f \n Total energy %f \n",
                    dt*(n+1),dt*timesteps,ek,ep,et);
     }
   }
  stop=time(NULL);
  printf("WallClock time %ld/n",stop-start)
   return 0;
}
