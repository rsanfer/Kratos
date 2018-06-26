//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes 

// External includes 

// Project includes
#include "custom_python/add_custom_processes_to_python.h"


// Processes
#include "custom_processes/adaptive_time_interval_process.hpp"
#include "custom_processes/split_elements_process.hpp"
#include "custom_processes/set_active_flag_process.hpp"

// Mesher initialization and finalization processes
#include "custom_processes/settle_fluid_model_structure_process.hpp"

// PreMeshing processes
#include "custom_processes/set_active_entities_mesher_process.hpp"
#include "custom_processes/recover_volume_losses_mesher_process.hpp"
#include "custom_processes/inlet_management_mesher_process.hpp"
#include "custom_processes/insert_new_nodes_mesher_process.hpp"
#include "custom_processes/remove_fluid_nodes_mesher_process.hpp"

// MiddleMeshing processes

// PostMeshing processes
#include "custom_processes/select_fluid_elements_mesher_process.hpp"


namespace Kratos
{
	
namespace Python
{
  	
void  AddCustomProcessesToPython(pybind11::module& m)
{

  using namespace pybind11;

  //**********MODEL STRUCTURE*********//
  class_<SettleFluidModelStructureProcess, SettleFluidModelStructureProcess::Pointer, SettleModelStructureProcess>
      (m, "FluidModelStructure")
      .def(init<ModelPart&, Flags, int>());

  
  //**********MESHER PROCESSES*********//
  
  class_<RemoveFluidNodesMesherProcess, RemoveFluidNodesMesherProcess::Pointer, MesherProcess>
      (m, "RemoveFluidNodes")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());
  
  class_<InsertNewNodesMesherProcess, InsertNewNodesMesherProcess::Pointer, MesherProcess>
      (m, "InsertNewNodes")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());
  
  class_<SelectFluidElementsMesherProcess, SelectFluidElementsMesherProcess::Pointer, MesherProcess>
      (m, "SelectFluidElements")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());

  class_<SetActiveEntitiesMesherProcess, SetActiveEntitiesMesherProcess::Pointer, MesherProcess>
      (m, "SetActiveEntities")
      .def(init<ModelPart&, bool, bool, int>());

  
  //*********SET SOLVER PROCESSES*************//

  class_<SetActiveFlagProcess, SetActiveFlagProcess::Pointer, MesherProcess>
      (m, "SetActiveFlagProcess")
      .def(init<ModelPart&, bool, bool, int>());

  class_<SplitElementsProcess, SplitElementsProcess::Pointer, Process>
      (m,"SplitElementsProcess")
      .def(init<ModelPart&, int>());
  

  //*********ADAPTIVE TIME STEP*************//
  
  class_<AdaptiveTimeIntervalProcess, AdaptiveTimeIntervalProcess::Pointer, Process>
      (m, "AdaptiveTimeIntervalProcess")
      .def(init<ModelPart&, int>());
  
  //***************INLET PROCESS************//
  
  class_<InletManagementProcess, InletManagementProcess::Pointer, Process>
      (m, "InletManagement")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());
  
  
  //*********VOLUME RECOVETY PROCESS********//
  
  class_<RecoverVolumeLossesMesherProcess, RecoverVolumeLossesMesherProcess::Pointer, MesherProcess>
      (m, "RecoverVolumeLosses")
      .def(init<ModelPart&, MesherUtilities::MeshingParameters&, int>());
  
}
 
}  // namespace Python.

} // Namespace Kratos

