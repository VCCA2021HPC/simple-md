/* A serial example molecular dynamics program
*
* compile using
* cc -lm md.c -o md
*
* run using
* ./md
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void forceupdate(float *position, float *force, 
                int numparticles, float sigma, float epsilon)
{

  // initialize to zero
  for(int i=0; i<3*numparticles;i++) force[i]=0;
  // add force components from other particles
  for(int i = 0; i < numparticles; i++)
  {
   for(int j = 0; j < i; j++)
    {
       float rsq = (position[3*i  ]-position[3*j  ])*(position[3*i  ]-position[3*j  ])
                  +(position[3*i+1]-position[3*j+1])*(position[3*i+1]-position[3*j+1])
                  +(position[3*i+2]-position[3*j+2])*(position[3*i+2]-position[3*j+2]);

      // loop over x,y and z components
      for(int k = 0; k < 3; k++)
      {
        force[3*i+k] +=
             24 * epsilon * ( (position[3*i+k]-position[3*j+k]) / (sigma*sigma) )  *
            ( 2 * pow(( sigma*sigma / rsq ),7) -  pow(( sigma*sigma / rsq ),4) );
      }
    }
   for(int j = i+1; j < numparticles; j++)
    {
       float rsq = (position[3*i  ]-position[3*j  ])*(position[3*i  ]-position[3*j  ])
                  +(position[3*i+1]-position[3*j+1])*(position[3*i+1]-position[3*j+1])
                  +(position[3*i+2]-position[3*j+2])*(position[3*i+2]-position[3*j+2]);

      // loop over x,y and z components
      for(int k=0;k<3;k++)
      {
        force[3*i+k] +=  
             24 * epsilon * ( (position[3*i+k] - position[3*j+k]) / (sigma*sigma) )  * 
            ( 2 * pow(( sigma*sigma / rsq ),7) -  pow(( sigma*sigma / rsq ),4) );
      }
    }
  }
  return;
}

void positionupdate(float *position, float *velocity, 
                    float dt, int numparticles)
{
    for(int i = 0; i<numparticles; i++)
    {
      position[i] +=  dt * velocity[i];
    }
    return;
}

void velocityupdate(float *velocity, float *force,
                    float dt, float mass, int numparticles)
{
    for(int i = 0; i < 3*numparticles; i++)
    {
         velocity[i] += dt * force[i] / mass;
    }
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
   return 0;
}
