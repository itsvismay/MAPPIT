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
			Matrix<double, 9, 9> Hi;

			for(int j=1; j<Q_i.rows() -1; j++){
				//for each node thats not the first node (0) and last node
				double q1 = Q_i.row(j-1)[0];
				double q2 = Q_i.row(j-1)[1];
				double q3 = Q_i.row(j-1)[2];
				double q4 = Q_i.row(j)[0];
				double q5 = Q_i.row(j)[1];
				double q6 = Q_i.row(j)[2];
				double q7 = Q_i.row(j+1)[0];
				double q8 = Q_i.row(j+1)[1];
				double q9 = Q_i.row(j+1)[2];
				double K = K_acc;
				//
				Vector3d x1 = Q_i.row(j-1);
				Vector3d v1 = Q_i.row(j) - Q_i.row(j-1);
				Vector3d v2 = Q_i.row(j+1) - Q_i.row(j);
				Vector3d v1xv2 = (v1.cross(v2)).transpose();
				double v1dv2 = v1.dot(v2);
				double v1xv2norm = v1xv2.norm();

				////----Dave's Energy-----
				double eps = 1e-5;
		        if(v1dv2<=eps){
		          v1dv2 = eps;
		          v2(0) += eps;
		          v2(1) += eps;
		          v2(2) += eps;
		        }
		       
		        v1xv2 = (v1.cross(v2)).transpose();
		        v1dv2 = v1.dot(v2);
		        v1xv2norm = v1xv2.norm();

		        double X = v1dv2;
		        double Y = v1xv2norm;

		        double left = 0;
		        double middle = 1.0/(1.0 + (Y/X)*(Y/X));
		        if(v1xv2.norm()>eps)
		          left = atan2(Y, X);

				//------------------Analytical gradient
		        Matrix3d v2CrossMatrix;
		        v2CrossMatrix<< 0,  v2(2), -v2(1),
		                       -v2(2),   0,  v2(0),
		                        v2(1), -v2(0),   0; 
		        Matrix3d v1CrossMatrix;
		        v1CrossMatrix<< 0,  -v1(2), v1(1),
		                         v1(2),   0,  -v1(0),
		                        -v1(1), v1(0),   0; 
				
				Vector3d right1;
				Vector3d right2;
				{
			        Vector3d lodhi1 = v1dv2*(1.0/v1xv2norm) *(v1xv2.transpose()*v2CrossMatrix);
			        Vector3d hidlo1 = v2*v1xv2norm;
			        Vector3d lodhi2 = v1dv2*(1.0/v1xv2norm) *(v1xv2.transpose()*v1CrossMatrix);
			        Vector3d hidlo2 = v1*v1xv2norm;
			        double lolo = v1dv2*v1dv2;
			        right1 = (lodhi1 - hidlo1)/lolo;
			        right2 = (lodhi2 - hidlo2)/lolo;
			    }
		        Vector3d g1 = left*middle*right1;
		        Vector3d g2 = left*middle*right2;


		        //--------------All the derivatives
		        //dleftdv1
		        Vector3d dleftdv1 = -middle*right1;

		        //dleftdv2
				Vector3d dleftdv2 = -middle*right2;

				//dmiddv1
				Vector3d dmiddv1 = (2*(Y/X))/(((Y/X)*(Y/X) + 1)*((Y/X)*(Y/X) + 1))*right1;

				//dmiddv2
				Vector3d dmiddv2 = (2*(Y/X))/(((Y/X)*(Y/X) + 1)*((Y/X)*(Y/X) + 1))*right2;

				//dright1dv1
				Matrix3d dright1dv1;
				{	
					Vector3d hi = (v1dv2*v1xv2.transpose()*v2CrossMatrix - v2.transpose()*v1xv2norm*v1xv2norm);
					double lo = v1dv2*v1dv2*v1xv2norm;
					Vector3d temp = (v1xv2.transpose()*v2CrossMatrix);
					// std::cout<<"v2m"<<std::endl;
					// std::cout<<v2CrossMatrix<<std::endl;
					// //v2m*v1xv2'*v2 - v1dv2*v2m'*v2m
					Matrix3d dhi_part1 = ((v2CrossMatrix*v1xv2)*v2.transpose() - v1dv2*v2CrossMatrix.transpose()*v2CrossMatrix);
					Matrix3d dhi_part2 = (v2*2*v1xv2.transpose()*v2CrossMatrix);
					Matrix3d dhi = dhi_part1 + dhi_part2;
					Matrix3d lodhi = lo*dhi;
					std::cout<<"hi"<<std::endl;
					std::cout<<hi<<std::endl;
					std::cout<<"lo"<<std::endl;
					std::cout<<lo<<std::endl;
					// std::cout<<"dhi part 1"<<std::endl;
					// std::cout<<dhi_part1<<std::endl;
					// std::cout<<"dhi part 2"<<std::endl;
					// std::cout<<dhi_part2<<std::endl;
					// std::cout<<"dhi"<<std::endl;
					// std::cout<<dhi<<std::endl;
					Vector3d dlo_part1 = -(2*v1dv2*v1xv2norm*v2);
					Vector3d dlo_part2 = -(v1dv2*v1dv2*(1.0/v1xv2norm) * (v1xv2.transpose()*v2CrossMatrix));
					Vector3d dlo = dlo_part1 + dlo_part2;
					std::cout<<"dlo part 1"<<std::endl;
					std::cout<<dlo_part1.transpose()<<std::endl;
					std::cout<<"dlo part 2"<<std::endl;
					std::cout<<dlo_part2.transpose()<<std::endl;
					std::cout<<"dlo"<<std::endl;
					std::cout<<dlo.transpose()<<std::endl;
					Matrix3d hidlo = hi*dlo.transpose();
					double lolo = lo*lo;
					dright1dv1 = (lodhi - hidlo)/lolo;
				}



				//dright2dv2
				Matrix3d dright2dv2;
				{	
					Vector3d hi = (v1dv2*v1xv2.transpose()*v1CrossMatrix - v1.transpose()*v1xv2norm*v1xv2norm);
					double lo = v1dv2*v1dv2*v1xv2norm;
					Matrix3d lodhi_part1 = (lo)*(v1*v1xv2.transpose()*v1CrossMatrix + v1dv2*v1CrossMatrix*v1CrossMatrix);
					Matrix3d lodhi_part2 = -(lo)*(v1*2*v1xv2.transpose()*v1CrossMatrix);
					Matrix3d lodhi = lodhi_part1 + lodhi_part2;
					Vector3d dlo_part1 = (2*v1dv2*v1xv2norm*v1);
					Vector3d dlo_part2 = (v1dv2*v1dv2*(1.0/v1xv2norm) * (v1xv2.transpose()*v1CrossMatrix));
					Vector3d dlo = dlo_part1 + dlo_part2;
					Matrix3d hidlo = hi*dlo.transpose();
					double lolo = lo*lo;
					dright2dv2 = (lodhi - hidlo)/lolo;
				}

		        //--------------

				
		      	Matrix3d H11;H11.setZero();// dGjp/dq(1:3)
		      	Matrix3d p1 = middle * dleftdv1 * right1.transpose(); //dleftdv*middle*right
		      	Matrix3d p2 = left * dmiddv1 * right1.transpose(); //left*dmiddledv1*right
				Matrix3d p3 = left * middle * dright1dv1; //left*middle*drightdv1
				H11 = -p1 + -p2 + -p3;
				
				std::cout<<"Left"<<std::endl;
				std::cout<<left<<std::endl;
				std::cout<<"Mid"<<std::endl;
				std::cout<<middle<<std::endl;
				std::cout<<"right1"<<std::endl;
				std::cout<<right1<<std::endl;
				std::cout<<"dleftdv1"<<std::endl;
				std::cout<<dleftdv1.transpose()<<std::endl;
				std::cout<<"dmiddv1"<<std::endl;
				std::cout<<dmiddv1.transpose()<<std::endl;
				std::cout<<"dright1dv1"<<std::endl;
				std::cout<<dright1dv1<<std::endl;
				std::cout<<"g1"<<std::endl;
				std::cout<<g1.transpose()<<std::endl;
				std::cout<<"##--P1"<<std::endl;
				std::cout<<p1<<std::endl;
				std::cout<<"##--P2"<<std::endl;
				std::cout<<p2<<std::endl;
				std::cout<<"##--P3"<<std::endl;
				std::cout<<p3<<std::endl;

				Matrix3d H12;H12.setZero();// dGjp/dq(4:6)
				
				Matrix3d H13;H13.setZero();// dGjp/dq(7:9)
				
				Matrix3d H31;H31.setZero();// dGjp/dq(1:3)
				H31 = H13.transpose();
				
				Matrix3d H32;H32.setZero();// dGjn/dq(1:3)
				
				Matrix3d H33;H33.setZero();// dGjp/dq(1:3)
				
				Matrix3d H21;H21.setZero();// dGjp/dq(1:3)
				H21 = H12.transpose();

				Matrix3d H22;H22.setZero();// dGjp/dq(1:3)

				Matrix3d H23;H23.setZero();// dGjp/dq(1:3)
				H23 = H32.transpose();

				Hi<<H11,H12, H13, H21, H22, H23, H31, H32, H33;

				// //---Etienne's Energy-------------------

				// double eps = 1e-5;
				//       if(v1xv2.norm()<=eps){
				//         v1xv2norm = eps;
				//       }
				// Vector3d z = v1xv2/v1xv2norm;
				// double X = v1.norm()*v2.norm() + v1dv2;
				// double Y = v1xv2.dot(z);
				// double angle = 2*atan2(Y, X);
				// double K = K_acc;

				// Vector3d Gjprev = -K * (angle) * (v1.cross(z)/v1.squaredNorm());
				// Vector3d Gjnext = -K * (angle) * (v2.cross(z)/v2.squaredNorm());
				// double z1 = z(0);
				// double z2 = z(1);
				// double z3 = z(2);

				// Matrix3d H11;H11.setZero();// dGjp/dq(1:3)
				
				// Matrix3d H12;H12.setZero();// dGjp/dq(4:6)
				
				// Matrix3d H13;H13.setZero();// dGjp/dq(7:9)
				
				// Matrix3d H31;H31.setZero();// dGjp/dq(1:3)
				// H31 = H13.transpose();
				
				// Matrix3d H32;H32.setZero();// dGjn/dq(1:3)
				
				// Matrix3d H33;H33.setZero();// dGjp/dq(1:3)
				
				// Matrix3d H21;H21.setZero();// dGjp/dq(1:3)
				// H21 = H12.transpose();

				// Matrix3d H22;H22.setZero();// dGjp/dq(1:3)
				
				// Matrix3d H23;H23.setZero();// dGjp/dq(1:3)
				// H23 = H32.transpose();

				// Hi<<H11,H12, H13, H21, H22, H23, H31, H32, H33;
				
				Vector3i inds;
				inds<<i*3*num_points_per_agent + 3*(j-1)+0, 
					i*3*num_points_per_agent + 3*j+0,
					i*3*num_points_per_agent + 3*(j+1)+0;

				for(int ii=0; ii<3; ii++)
				{
					for(int jj=0; jj<3; jj++)
					{
						for(int kk=0; kk<3; kk++)
						{
							for(int ll=0; ll<3; ll++)
							{
								H(inds(ii)+kk, inds(jj)+ll) += Hi(3*ii+kk, 3*jj+ll);
							}
						}
						
					}
				}


				//std::cout<<"------"<<std::endl;
				//std::cout<<Hi<<std::endl;

				//std::cout<<"------ Check symmetric"<<std::endl;
				//MatrixXd Hit = Hi.transpose();
				//std::cout<< Hi - Hit<<std::endl;
			}
		}

	}
}