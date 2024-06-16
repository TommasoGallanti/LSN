#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include <cassert>
using namespace std;


class Posizione {
   
   public:
        
    Posizione() : m_pos(3), m_cont(3){
            
    }

     Posizione(double x, double y, double z) :m_pos(3) , m_cont(3){
         m_pos[0]=x ;
         m_pos[1]=y ;
         m_pos[2]=z ;
     }

     //distruttore
     ~Posizione() {}

     //metodi
     double getX() const { return m_pos[0]; }       //coord cartesiane
     double getY() const { return m_pos[1]; }
     double getZ() const { return m_pos[2]; } 
     double getR() const {return sqrt(m_pos[0]*m_pos[0]+m_pos[1]*m_pos[1]+m_pos[2]*m_pos[2]); }   //coord sferiche 
     double getPhi() const {return atan2(m_pos[1],m_pos[0]); }
     double getTheta() const {return acos(m_pos[2]/getR()); }
     double getRho() const {return sqrt(m_pos[0]*m_pos[0]+m_pos[1]*m_pos[1]); }          //raggio coord cilindriche 
     double Distanza(const Posizione& b) const {                       //distanza punti
         return sqrt ( pow(getX()-b.getX(),2) +
                       pow(getY()-b.getY(),2) +
                       pow(getZ()-b.getZ(),2) ); 
     } 
     void SetStepLenght(double a) { m_a = a;}
     double GetPos(int i) const { 
            assert(i>=0 && i<3);
            return m_pos[i];}

     void SetPos(int i, double x) {
            assert(i>=0 && i<3) ;
            m_pos[i] = x;}

     void AddPos(int i, double x) {
            assert(i>=0 && i<3) ;
            m_pos[i] += x;}

     void CopyPos( Posizione p) { 
            for(int i{}; i<3; ++i) 
                m_pos[i] = p.GetPos(i) ;
            }
    
    

     double GetStepLenght() {return m_a;}
     double GetR2(){ return (m_pos[0]*m_pos[0] +m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]);}
     void SetR2(){m_r2 =(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]); }
        
    
   protected:
        
     vector<double> m_pos;   // pos lattice x, y, z
     vector<double> m_cont;  //pos continue, r,theta,phi
     double m_a = 1.;        //passo
     double m_r2;            //raggio^2 cont

};

class Walk : public Posizione {

    public:

        Walk() : m_step {100} , traj(m_step) {}
        Walk(int step) : m_step {step} , traj(m_step) {}   
        

        void Random_Walk(vector<double> rnd) ;   //generate a discrete random walk 
        void Random_Cont(vector<double> rnd);     //generate cont rw 
        

        void Settraj(int i, Posizione p) { 
            assert(i>=0 && i< int( traj.size() ) );
            traj[i].CopyPos(p);
        }

        double GetR2(int i) {
            return traj[i].GetR2();
        }

        void SetR2(int i) {
            traj[i].SetR2();
        }

        double GetStep() {return m_step;}

        void Print() {
            for(Posizione value : traj) {
                for(int i{}; i< 3; ++i)
                    cout<<value.GetPos(i)<<"\t"<<endl;
            }
        }
       
        Walk& operator= ( const Walk& w) {
            m_step = w.m_step;
            traj.resize(0); 

            for(Posizione value : w.traj) {
                traj.push_back(value);
            }
            return *this;
        }
        
      
    private:

        int m_step;
        vector<Posizione> traj;

};

void Walk::Random_Walk(vector<double> rnd) {     //discrete random Walk 

    Posizione v1, v2;

    for(size_t i{}; i< traj.size(); ++i) {
        double rnd_step = 3*rnd[i]; 
        double rnd_sign = 1 + rnd[traj.size() + i];
        double sign = (rnd_sign > 1.5 ? 1 : -1);

        for(int j{}; j < 3; ++j) {
            v1.SetPos(j, 0.);
        }

        v1.SetPos(int(rnd_step), sign*v1.GetStepLenght());

        for(size_t j{}; j < m_pos.size(); ++j) {
            v2.AddPos(j, v1.GetPos(j));
        }

        Settraj(i, v2);
        GetR2(i);       //charge the value of r2 in each Walk Posizione
    }
}

        
void Walk::Random_Cont(vector<double> rnd) {   //cont RW
            
        Posizione v1, v2;

        for(size_t i{}; i< traj.size(); ++i) {

                
            for(int i{}; i < 3; ++i)
                        v1.SetPos(i,0.);

            double rnd_theta = 2*M_PI*rnd[i]; 
            double rnd_phi = acos( 1 -2*rnd[traj.size() + i] ); //Sampling of a solid angle
                
            v1.SetPos(0 , m_a*sin(rnd_phi)*cos(rnd_theta)  ); //x
            v1.SetPos(1 , m_a*sin(rnd_phi)*sin(rnd_theta)  ); //y
            v1.SetPos(2 , m_a*cos(rnd_phi)  );                //z

            for(size_t i{}; i < m_pos.size(); ++i )      
            v2.AddPos(i,v1.GetPos(i));


            Settraj(i,v2);
            SetR2(i); //charge the value of  r2 in each Walk Posizione
        }
} 






















