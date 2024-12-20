// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Yixuan Tang
// =============================================================================
//
// FEA for 3D beams of 'cable' type (ANCF gradient-deficient beams)
//       uses the Chrono MKL module
//
// =============================================================================

#include <iomanip>
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/physics/ChLinkMotorAll.h"

#include "chrono/timestepper/ChTimestepper.h"
#include "chrono_pardisomkl/ChSolverPardisoMKL.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"

#include "FEAcables.h"
///////////////
#include "chrono/fea/ChElementCableANCF.h"
#include "chrono/fea/ChElementGeneric.h"
#include "chrono/functions/ChFunctionSine.h"
#include "chrono/fea/ChMesh.h"




using namespace chrono;
using namespace fea;

const std::string out_dir = GetChronoOutputPath() + "ANCFCable";  //----D:\Chrono\chrono_build\bin

bool cantilever = false;
bool is_already_disabled = false;
    // ChSystemSMC sys;
 
int main(int argc, char* argv[]) {
    std::cout << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << std::endl;

    // Create a Chrono::Engine physical system
    ChSystemSMC sys;
  
    sys.SetNumThreads(std::min(4, ChOMP::GetNumProcs()), 0, 1);
    //sys.SetNumThreads(1, 0, 1);

    // Create a mesh, that is a container for groups of elements and
    // their referenced nodes.
    auto my_mesh = chrono_types::make_shared<ChMesh>();
   

   auto msection_cable2 = chrono_types::make_shared<ChBeamSectionCable>();

   //// This ChBuilderCableANCF helper object is very useful because it will
   //// subdivide 'beams' into sequences of finite elements of beam type, ex.
   //// one 'beam' could be made of 5 FEM elements of ChElementBeamANCF_3333 class.
   //// If new nodes are needed, it will create them for you.
   ChBuilderCableANCF builder;
   
   //// Now, simply use BuildBeam to create a beam from a point to another:
   if (!cantilever) {
       msection_cable2->SetDiameter(0.4);
       msection_cable2->SetYoungModulus(10.88e10);
       msection_cable2->SetRayleighDamping(0.1);
       builder.BuildBeam(my_mesh,          // the mesh where to put the created nodes and elements
                         msection_cable2,  // the ChBeamSectionCable to use for the ChElementBeamANCF_3333 elements
                         50,               // the number of ChElementBeamANCF_3333 to create
                         ChVector3d(-837.8, -200, 0),   // the 'A' point in space (beginning of beam)  837.8, -200
                         ChVector3d(0, 0, 0));  // the 'B' point in space (end of beam)
       my_mesh->SetAutomaticGravity(true);
         }
 else{
       msection_cable2->SetDiameter(0.01);
       msection_cable2->SetYoungModulus(10.88e10);
       msection_cable2->SetRayleighDamping(0.1);
       builder.BuildBeam(my_mesh,          // the mesh where to put the created nodes and elements
                         msection_cable2,  // the ChBeamSectionCable to use for the ChElementBeamANCF_3333 elements
                         10,               // the number of ChElementBeamANCF_3333 to create
                         ChVector3d(0, 0, 0),  // the 'A' point in space (beginning of beam)  837.8, -200
                         ChVector3d(2, 0, 0));     // the 'B' point in space (end of beam)
       my_mesh->SetAutomaticGravity(false);
   }

   //// For instance, now retrieve the A end and add a constraint to
   //// block the position only of that node:
      
   //ChVector3d dir1(1, 0, 0);
   //auto ground = chrono_types::make_shared<ChNodeFEAxyzD>(ChVector3d(0, 0, 0), dir1);
   //ground->SetFixed(true);
   //my_mesh->AddNode(ground);


   auto ground = chrono_types::make_shared<ChBody>();
   ground->SetFixed(true);
   sys.Add(ground);


   auto moving_body = chrono_types::make_shared<ChBody>();
   auto constraint_prismatic = chrono_types::make_shared<ChLinkMotorLinearPosition>();
   auto constraint_hinge1 = chrono_types::make_shared<ChLinkNodeFrame>();
   auto ctrl_fun = chrono_types::make_shared<ChFunctionSine>();

   if (!cantilever) {
       sys.Add(moving_body);    
       constraint_prismatic->Initialize(moving_body, ground,
                                       ChFramed(builder.GetLastBeamNodes().back()->GetPos(), Q_ROTATE_Z_TO_Y));
       sys.Add(constraint_prismatic);
       constraint_hinge1->Initialize(builder.GetLastBeamNodes().back(), moving_body);
       sys.Add(constraint_hinge1);
   }

   std::shared_ptr<ChNodeFEAxyzD> node_A;
    node_A = builder.GetLastBeamNodes().front();
   auto constraintxyz = chrono_types::make_shared<ChLinkNodeFrame>();
   constraintxyz->Initialize(node_A, ground);
   sys.Add(constraintxyz);
   node_A->SetSlope1Fixed(true);
   double x = sys.GetChTime();
   
    // Remember to add the mesh to the system!
    sys.Add(my_mesh);

    //// Visualization of the FEM mesh.
    //// This will automatically update a triangle mesh (a ChVisualShapeTriangleMesh
    //// asset that is internally managed) by setting  proper
    //// coordinates and vertex colors as in the FEM elements.
    //// Such triangle mesh can be rendered by Irrlicht or POVray or whatever
    //// postprocessor that can handle a colored ChVisualShapeTriangleMesh).

    //auto mvisualizebeamA = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    //mvisualizebeamA->SetFEMdataType(ChVisualShapeFEA::DataType::ELEM_BEAM_MZ);
    //mvisualizebeamA->SetColorscaleMinMax(-0.4, 0.4);
    //mvisualizebeamA->SetSmoothFaces(true);
    //mvisualizebeamA->SetWireframe(false);
    //my_mesh->AddVisualShapeFEA(mvisualizebeamA);

    //auto mvisualizebeamC = chrono_types::make_shared<ChVisualShapeFEA>(my_mesh);
    //mvisualizebeamC->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_CSYS);
    //mvisualizebeamC->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    //mvisualizebeamC->SetSymbolsThickness(0.006);
    //mvisualizebeamC->SetSymbolsScale(0.01);
    //mvisualizebeamC->SetZbufferHide(false);
    //my_mesh->AddVisualShapeFEA(mvisualizebeamC);

    //// Create the Irrlicht visualization system
    //auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    //vis->SetWindowSize(800, 600);
    //vis->SetWindowTitle("Cables FEM (MKL)");
    //vis->Initialize();
    //vis->AddLogo();
    //vis->AddSkyBox();
    //vis->AddTypicalLights();

    //if (!cantilever) {
    //    vis->AddCamera(ChVector3d(400, 0, 200.0));
    //}
    //else vis->AddCamera(ChVector3d(0, 0, 2.5));
    //vis->AttachSystem(&sys);

    // Configure PardisoMKL solver.
    // For this simple and relatively small problem, use of the sparsity pattern learner may introduce additional
    // overhead (if the sparsity pattern is not locked).

    auto mkl_solver = chrono_types::make_shared<ChSolverPardisoMKL>(1);
    mkl_solver->UseSparsityPatternLearner(false);
    mkl_solver->LockSparsityPattern(false);
    mkl_solver->SetVerbose(false);
    sys.SetSolver(mkl_solver);

    sys.Update();

    // Set integrator
    sys.SetTimestepperType(ChTimestepper::Type::HHT);
   //  sys.SetTimestepperType(ChTimestepper::Type::EULER_IMPLICIT_LINEARIZED);
    auto ts = chrono_types::make_shared<ChTimestepperHHT>(&sys);
    sys.SetTimestepper(ts);
    ts->SetStepControl(false);

    double N_ELE_TEST0 = 9;  // the last element is the ploating point
    double N_ELE_TEST1 = 24;
    double N_ELE_TEST2 = 39;


    auto element0 = builder.GetLastBeamElements().at(N_ELE_TEST0); //---choosed
    std::shared_ptr<ChNodeFEAxyzD> node_element0;
    node_element0 = builder.GetLastBeamNodes().at(N_ELE_TEST0);//  at(9) Origin 

    auto element1 = builder.GetLastBeamElements().at(N_ELE_TEST1);  //---choosed
    std::shared_ptr<ChNodeFEAxyzD> node_element1;
    node_element1 = builder.GetLastBeamNodes().at(N_ELE_TEST1);  //  at(24) Origin 

    auto element2 = builder.GetLastBeamElements().at(N_ELE_TEST2);  //---choosed
    std::shared_ptr<ChNodeFEAxyzD> node_element2;
    node_element2 = builder.GetLastBeamNodes().at(N_ELE_TEST2);  //  at(39) Origin 

    std::shared_ptr<ChNodeFEAxyzD> node_floating;
    node_floating = builder.GetLastBeamNodes().back();
    // node_floating = builder.GetLastBeamNodes().front();


    sys.DoAssembly(0);
    constexpr double time_step = 0.01;
    constexpr double time_step_prt = 0.01;
    constexpr double time_length =100;
    
    int Nframes = (int)(time_length / time_step);
    int itv_frame = (int)(time_step_prt / time_step);

    int frame = 0;
    int Nframe = 0;

   //ChVectorDynamic<> Fi0(12);
   //ChMatrixDynamic<> ForceGlobal;
   //ForceGlobal.resize(Nframes, 6);
   //ForceGlobal.setZero();

    int STEP_SIZE = 10;

   ChMatrixDynamic<> Pos;
   Pos.resize(Nframes/STEP_SIZE , 27); // 9*3
   Pos.setZero();

   ChVector3d Fforce0;
   ChVector3d Mmoment0;
   ChVector3d Fforce1;
   ChVector3d Mmoment1;
   ChVector3d Fforce2;
   ChVector3d Mmoment2;
   ChMatrixDynamic<> ForceMomentSec;
   ForceMomentSec.resize(Nframes / STEP_SIZE, 18);  // 6*3
   ForceMomentSec.setZero();



 //  ChVectorDynamic<> Fi(sys.GetNumCoordsPosLevel()); X FOR WHOLE SYSTEM

  // msection_cable2;
   // m_section;
    // Simulation loop
    ChTimer m_timer_computation;
    double time = 0;
     
    //  if (!cantilever) {
    //    auto ctrl_fun = chrono_types::make_shared<ChFunctionSine>();
    //    ctrl_fun->SetAmplitude(2);  // set amplitude;
    //    ctrl_fun->SetFrequency(1);  // set frequency; T=1;
    //    constraint_prismatic->SetMotionFunction(ctrl_fun);
    //} 
    //  else {
    //    auto ctrl_fun = chrono_types::make_shared<ChFunctionSine>();
    //    ctrl_fun->SetAmplitude(0.1);  // set amplitude;
    //    ctrl_fun->SetFrequency(1);    // set frequency; T=1;
    //    constraint_prismatic->SetMotionFunction(ctrl_fun);
    //}

    // while (vis->Run()) {

    while (frame < Nframes) {
        //m_timer_computation.start();
        //vis->BeginScene();
        //vis->Render();

        if (cantilever) {
            if (sys.GetChTime() < 0.2) {
                node_floating->SetForce(ChVector3d(0, -0.2, 0));
            } else {
                node_floating->SetForce(ChVector3d(0, 0, 0));
            }
        } else {
            if (sys.GetChTime() > 0.5) {
                ctrl_fun->SetAmplitude(0);  // set amplitude;
                ctrl_fun->SetFrequency(1);  // set frequency; T=1;
                constraint_prismatic->SetMotionFunction(ctrl_fun);
            } else {
                ctrl_fun->SetAmplitude(2);  // set amplitude;
                ctrl_fun->SetFrequency(1);  // set frequency; T=1;
                constraint_prismatic->SetMotionFunction(ctrl_fun);
            }
        }

        sys.DoStepDynamics(time_step);
        //vis->EndScene();

        //m_timer_computation.stop();
        //time= m_timer_computation.GetTimeSeconds();
        //m_timer_computation.reset();

        //model.PrintBodyPositions();
       // ChVector3d force = builder.GetLastBeamNodes()[0]->GetForce();
        if (frame % (itv_frame * STEP_SIZE) == 0) {
            // element0->ComputeInternalForces(Fi0); // This is the global force
            element0->EvaluateSectionForceTorque(0, Fforce0, Mmoment0); 

            element1->EvaluateSectionForceTorque(0, Fforce1, Mmoment1); 

            element2->EvaluateSectionForceTorque(0, Fforce2, Mmoment2); 

            if (frame % (itv_frame * STEP_SIZE) == 0) {
                std::cout << "Time" << frame * time_step << "\n";
             }
                  //ForceGlobal(Nframe, 0) = Fi0(0);
                  //ForceGlobal(Nframe, 1) = Fi0(1);
                  //ForceGlobal(Nframe, 2) = Fi0(2);
                  //ForceGlobal(Nframe, 3) = Fi0(3);
                  //ForceGlobal(Nframe, 4) = Fi0(4);
                  //ForceGlobal(Nframe, 5) = Fi0(5);
                  //-----------------------------------------------------------------------
                  Pos(Nframe, 0) = node_element0->GetPos().x();
                  Pos(Nframe, 1) = node_element0->GetPos().y();
                  Pos(Nframe, 2) = node_element0->GetPos().z();
                  Pos(Nframe, 3) = node_element0->GetPosDt().x();
                  Pos(Nframe, 4) = node_element0->GetPosDt().y();
                  Pos(Nframe, 5) = node_element0->GetPosDt().z();
                  Pos(Nframe, 6) = node_element0->GetPosDt2().x();
                  Pos(Nframe, 7) = node_element0->GetPosDt2().y();
                  Pos(Nframe, 8) = node_element0->GetPosDt2().z();

                  Pos(Nframe, 9) = node_element1->GetPos().x();
                  Pos(Nframe, 10) = node_element1->GetPos().y();
                  Pos(Nframe, 11) = node_element1->GetPos().z();
                  Pos(Nframe, 12) = node_element1->GetPosDt().x();
                  Pos(Nframe, 13) = node_element1->GetPosDt().y();
                  Pos(Nframe, 14) = node_element1->GetPosDt().z();
                  Pos(Nframe, 15) = node_element1->GetPosDt2().x();
                  Pos(Nframe, 16) = node_element1->GetPosDt2().y();
                  Pos(Nframe, 17) = node_element1->GetPosDt2().z();

                  Pos(Nframe, 18) = node_element2->GetPos().x();
                  Pos(Nframe, 19) = node_element2->GetPos().y();
                  Pos(Nframe, 20) = node_element2->GetPos().z();
                  Pos(Nframe, 21) = node_element2->GetPosDt().x();
                  Pos(Nframe, 22) = node_element2->GetPosDt().y();
                  Pos(Nframe, 23) = node_element2->GetPosDt().z();
                  Pos(Nframe, 24) = node_element2->GetPosDt2().x();
                  Pos(Nframe, 25) = node_element2->GetPosDt2().y();
                  Pos(Nframe, 26) = node_element2->GetPosDt2().z();

             /*   Pos(Nframe, 27) = node_floating->GetPos().x();
                  Pos(Nframe, 28) = node_floating->GetPos().y();
                  Pos(Nframe, 29) = node_floating->GetPos().z();*/
                  //-----------------------------------------------------------------------
                  ForceMomentSec(Nframe, 0) = Fforce0.x();
                  ForceMomentSec(Nframe, 1) = Fforce0.y();
                  ForceMomentSec(Nframe, 2) = Fforce0.z();
                  ForceMomentSec(Nframe, 3) = Mmoment0.x();
                  ForceMomentSec(Nframe, 4) = Mmoment0.y();
                  ForceMomentSec(Nframe, 5) = Mmoment0.z();

                  ForceMomentSec(Nframe, 6) = Fforce1.x();
                  ForceMomentSec(Nframe, 7) = Fforce1.y();
                  ForceMomentSec(Nframe, 8) = Fforce1.z();
                  ForceMomentSec(Nframe, 9) = Mmoment1.x();
                  ForceMomentSec(Nframe, 10) = Mmoment1.y();
                  ForceMomentSec(Nframe, 11) = Mmoment1.z();

                  ForceMomentSec(Nframe, 12) = Fforce2.x();
                  ForceMomentSec(Nframe, 13) = Fforce2.y();
                  ForceMomentSec(Nframe, 14) = Fforce2.z();
                  ForceMomentSec(Nframe, 15) = Mmoment2.x();
                  ForceMomentSec(Nframe, 16) = Mmoment2.y();
                  ForceMomentSec(Nframe, 17) = Mmoment2.z();

                  std::ofstream file_defl_corot1((out_dir + "/MooringLineD0.001(M)0.5_100s_h=0.005_Pos1.dat").c_str());
                  std::ofstream file_defl_corot2((out_dir + "/MooringLineD0.001(M)0.5_100s_h=0.005_ForceMomentSec1.dat").c_str());

                 file_defl_corot1 << std::setprecision(12) << std::scientific;
                 file_defl_corot1 << Pos;
                 file_defl_corot2 << std::setprecision(12) << std::scientific;
                 file_defl_corot2 << ForceMomentSec;

                 Nframe++;
           
        }
        frame++;

    }

    return 0;
}
