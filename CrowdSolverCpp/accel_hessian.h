#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace Eigen;

namespace crowds{
	void accel_hessian(VectorXd& q, int num_agents, int num_points_per_agent, double K_acc, MatrixXd& H)
	{
		H.resize(q.size(), q.size());
		H.setZero();

		double e = 0.0;
		for(int i=0; i<num_agents; i++)
		{
			double mass_i = 1.0;//hard coded for now

			VectorXd q_i = q.segment(i*3*num_points_per_agent, 3*num_points_per_agent);
			MatrixXd Q_i = Map<MatrixXd>(q_i.data(), 3, q_i.size()/3).transpose();


			for(int j=1; j<Q_i.rows() -1; j++){
				//for each node thats not the first node (0) and last node
				Vector3d x1 = Q_i.row(j-1);
				Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
				Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
				Vector3d v1xv2 = (v1.cross(v2)).transpose();
				double v1dv2 = v1.dot(v2);
				double v1xv2norm = v1xv2.norm();
				double eps = 1e-3;
		        if(v1xv2.norm()<=eps){
		          v1xv2norm = eps;
		        }
				Vector3d z = v1xv2/v1xv2norm;
				double X = v1.norm()*v2.norm() + v1dv2;
				double Y = v1xv2.dot(z);
				double angle = 2*atan2(Y, X);
				double K = K_acc;

				Vector3d Gjprev = -K * (angle) * (v1.cross(z)/v1.squaredNorm());
        		Vector3d Gjnext = -K * (angle) * (v2.cross(z)/v2.squaredNorm());
        
        		Matrix3d cross_z;
        		cross_z<< 0,     z[2], -z[1],
						 -z[2],     0,  z[0],
						  z[1], -z[0],     0;
				cross_z*= -1;

				Vector3d v1xz = v1.cross(z);
				Vector3d v1xz_nv1 = (v1.cross(z)/v1.squaredNorm());
				Vector3d v2xz = v2.cross(z);
				Vector3d v2xz_nv2 = (v2.cross(z)/v2.squaredNorm());

        		//H11 -> d Gjprev/d xi-1
				Matrix3d H11 = K*v1xz_nv1*v1xz_nv1.transpose() - K*angle*((-v1.squaredNorm()*cross_z + 2*v1*v1xz.transpose())/(v1.squaredNorm()*v1.squaredNorm())).transpose();
				//H12 -> d Gjprev/d xi
				Matrix3d H12 = K*v1xz_nv1*(v1xz_nv1+v2xz_nv2).transpose() - K*angle*((v1.squaredNorm()*cross_z - 2*v1*v1xz.transpose())/(v1.squaredNorm()*v1.squaredNorm()));
				//H13 -> d Gjprev/d xi+1
				Matrix3d H13 = K*v1xz_nv1*v2xz_nv2.transpose();

				//H31 -> d Gjnext/d xi-1
				Matrix3d H31 = H13.transpose();
				//H32 -> d Gjnext/d xi
				Matrix3d H32; K*v2xz_nv2*(v1xz_nv1+v2xz_nv2).transpose() - K*angle*((-v2.squaredNorm()*cross_z + 2*v2*v2xz.transpose())/(v2.squaredNorm()*v2.squaredNorm()));
				//H33 -> d Gjnext/d xi+1
				Matrix3d H33 = K*v2xz_nv2*v2xz_nv2.transpose() - K*angle*((v2.squaredNorm()*cross_z  - 2*v2*v2xz.transpose())/(v2.squaredNorm()*v2.squaredNorm()));


				//H21 -> d (Gjprev + Gjnext)/d xi
				Matrix3d H21 = H12.transpose();

				//H21 -> d (Gjprev + Gjnext)/d xi
				Matrix3d H22 = H12 + H32;

				//H21 -> d (Gjprev + Gjnext)/d xi+1
				Matrix3d H23 = H13 + H33 ;


				

				std::cout<<"------"<<std::endl;
				H<<H11,H12, H13, H21, H22, H23, H31, H32, H33;
				std::cout<<H<<std::endl;

				std::cout<<"------ Check symmetric"<<std::endl;
				MatrixXd Ht = H.transpose();
				std::cout<< H - Ht<<std::endl;
			}
		}

	}
}