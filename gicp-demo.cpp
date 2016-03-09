// g2o - General Graph Optimization
// Copyright (C) 2011 Kurt Konolige
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
//   notice, this list of conditions and the following disclaimer in the
//   documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Eigen/StdVector>
#include <random>
#include <iostream>
#include <stdint.h>

#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/types/icp/types_icp.h"

using namespace Eigen;
using namespace std;
using namespace g2o;

// sampling distributions
  class Sample
  {

    static default_random_engine gen_real;
    static default_random_engine gen_int;
  public:
    static int uniform(int from, int to);

    static double uniform();

    static double gaussian(double sigma);
  };


  default_random_engine Sample::gen_real;
  default_random_engine Sample::gen_int;

  int Sample::uniform(int from, int to)
  {
    uniform_int_distribution<int> unif(from, to);
    int sam = unif(gen_int);
    return  sam;
  }

  double Sample::uniform()
  {
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    double sam = unif(gen_real);
    return  sam;
  }

  double Sample::gaussian(double sigma)
  {
    std::normal_distribution<double> gauss(0.0, sigma);
    double sam = gauss(gen_real);
    return  sam;
  }


Eigen::Affine3d toAffine(const Eigen::Isometry3d& pose)
{
	/*--------------------------------------------
	*			convert to Eigen::Affine3d
	*--------------------------------------------*/	
	Eigen::Affine3d p(pose.rotation());
	p.translation() = pose.translation();

	return p;
}



void addTransformationEdge(g2o::SparseOptimizer* graph, g2o::OptimizableGraph::Vertex *from, 
	g2o::OptimizableGraph::Vertex *to, 
	const Eigen::Isometry3d transform, 
	int edge_id )
{

	assert(from != 0 && to != 0);
//	ROS_INFO_STREAM("Adding adge from: "<<from->id()<<" to : "<< to->id());
//	g2o::RobustKernelTukey* tk = new g2o::RobustKernelTukey;
	g2o::EdgeSE3* edge = new g2o::EdgeSE3();
	edge->setId(edge_id);
	edge->resize(2);
	edge->setVertex(0, from);
	edge->setVertex(1, to);
	edge->setMeasurement(transform);
	//edge->setInformation(dvo::core::Matrix6d::Identity()*0.008*0.008);
	edge->setInformation(Eigen::Matrix<double,6,6>::Identity() );
//	edge->setRobustKernel(tk);
	graph->addEdge(edge);


}

//
// set up simulated system with noise, optimize it
//

int main()
{
  double euc_noise = 0.01;       // noise in position, m
  //  double outlier_ratio = 0.1;


  SparseOptimizer optimizer;
  optimizer.setVerbose(false);

  // variable-size block solver
  BlockSolverX::LinearSolverType * linearSolver = new LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
  BlockSolverX * solver_ptr = new BlockSolverX(linearSolver);
  g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);

  optimizer.setAlgorithm(solver);

  /*
  vector<Vector3d> true_points;
  for (size_t i=0;i<1000; ++i)
  {
    true_points.push_back(Vector3d((Sample::uniform()-0.5)*3,
                                   Sample::uniform()-0.5,
                                   Sample::uniform()+10));
  }*/


  // set up two poses
  int vertex_id = 0;
  for (size_t i=0; i<4; ++i)
  {
    // set up rotation and translation for this node
    Vector3d t(0,0,i);
    Quaterniond q;
    q.setIdentity();

    Eigen::Isometry3d cam; // camera pose
    cam = q;
    cam.translation() = t;

    // set up node
    VertexSE3 *vc = new VertexSE3();
    vc->setEstimate(cam);

    vc->setId(vertex_id);      // vertex id

    cerr << t.transpose() << " | " << q.coeffs().transpose() << endl;

    // set first cam pose fixed
    if (i==0)
      vc->setFixed(true);

    // add to optimizer
    optimizer.addVertex(vc);

    vertex_id++;                
  }
  // set up edges for these
  // create transform matrix 
  Vector3d t(0,0,2);
  t += Vector3d(Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ));
  Quaterniond q;
  q.setIdentity();
  Eigen::Isometry3d edge_cam;
  edge_cam = q;
  edge_cam.translation() = t;
  Eigen::Affine3d edge_tr = toAffine(edge_cam);

  addTransformationEdge(&optimizer,optimizer.vertex(0), optimizer.vertex(2), edge_cam, 3);

  Vector3d t2(0,0,1);
  t2 += Vector3d(Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ));
  Quaterniond q2;
  q2.setIdentity();
  Eigen::Isometry3d edge_cam2;
  edge_cam2 = q2;
  edge_cam2.translation() = t2;
  addTransformationEdge(&optimizer,optimizer.vertex(0), optimizer.vertex(1), edge_cam2, 1);

  Vector3d t3(0,0,1);
  t3 += Vector3d(Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ));
  Quaterniond q3;
  q3.setIdentity();
  Eigen::Isometry3d edge_cam3;
  edge_cam3 = q3;
  edge_cam3.translation() = t3;
  
  addTransformationEdge(&optimizer,optimizer.vertex(1), optimizer.vertex(2), edge_cam3, 2);


  /*
  // set up point matches
  for (size_t i=0; i<true_points.size(); ++i)
  {
    // get two poses
    VertexSE3* vp0 = 
      dynamic_cast<VertexSE3*>(optimizer.vertices().find(0)->second);
    VertexSE3* vp1 = 
      dynamic_cast<VertexSE3*>(optimizer.vertices().find(1)->second);

    // calculate the relative 3D position of the point
    Vector3d pt0,pt1;
    pt0 = vp0->estimate().inverse() * true_points[i];
    pt1 = vp1->estimate().inverse() * true_points[i];

    // add in noise
    pt0 += Vector3d(Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ));

    pt1 += Vector3d(Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ),
                    Sample::gaussian(euc_noise ));

    // form edge, with normals in varioius positions
    Vector3d nm0, nm1;
    nm0 << 0, i, 1;
    nm1 << 0, i, 1;
    nm0.normalize();
    nm1.normalize();

    Edge_V_V_GICP * e           // new edge with correct cohort for caching
        = new Edge_V_V_GICP(); 

    e->setVertex(0, vp0);      // first viewpoint

    e->setVertex(1, vp1);      // second viewpoint

    EdgeGICP meas;
    meas.pos0 = pt0;
    meas.pos1 = pt1;
    meas.normal0 = nm0;
    meas.normal1 = nm1;

    e->setMeasurement(meas);
    //        e->inverseMeasurement().pos() = -kp;
    
    meas = e->measurement();
    // use this for point-plane
    e->information() = meas.prec0(0.01);

    // use this for point-point 
    //    e->information().setIdentity();

    //    e->setRobustKernel(true);
    //e->setHuberWidth(0.01);

    optimizer.addEdge(e);
  }*/

  // move second cam off of its true position
  VertexSE3* vc = 
    dynamic_cast<VertexSE3*>(optimizer.vertices().find(1)->second);
  Eigen::Isometry3d cam = vc->estimate();
  cam.translation() = Vector3d(0,0,0.2);
  vc->setEstimate(cam);

  VertexSE3* vc2 = 
    dynamic_cast<VertexSE3*>(optimizer.vertices().find(2)->second);
  Eigen::Isometry3d cam2 = vc2->estimate();
  cam2.translation() = Vector3d(0,0,0.2);
  vc2->setEstimate(cam2);

  cout << endl << "Second vertex should be near 0,0,1" << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(0)->second)
    ->estimate().translation().transpose() << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(1)->second)
    ->estimate().translation().transpose() << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(2)->second)
    ->estimate().translation().transpose() << endl;

  optimizer.initializeOptimization();
  optimizer.computeActiveErrors();
  cout << "Initial chi2 = " << FIXED(optimizer.chi2()) << endl;

  optimizer.setVerbose(true);

  optimizer.optimize(50);

  cout << endl << "Second vertex should be near 0,0,1" << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(0)->second)
    ->estimate().translation().transpose() << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(1)->second)
    ->estimate().translation().transpose() << endl;
  cout <<  dynamic_cast<VertexSE3*>(optimizer.vertices().find(2)->second)
    ->estimate().translation().transpose() << endl;
}
