#include "flecsaleAdaptor.h"


namespace FlecsaleAdaptor
{

vtkCPProcessor* Processor = NULL;
vtkUnstructuredGrid* VTKGrid;


void Init(int numScripts, char* scripts[])
{
	if (Processor == NULL)
	{
		Processor = vtkCPProcessor::New();
		Processor->Initialize();
	}
	else
	{
		Processor->RemoveAllPipelines();
	}
  

	for (int i=1; i<numScripts; i++)
	{
		vtkNew<vtkCPPythonScriptPipeline> pipeline;
		pipeline->Initialize(scripts[i]);
		Processor->AddPipeline(pipeline.GetPointer());
	}
}


void Finalize()
{
	if (Processor)
	{
		Processor->Delete();
		Processor = NULL;
	}
	if (VTKGrid)
	{
		VTKGrid->Delete();
		VTKGrid = NULL;
	}
}



void CoProcess(vtkUnstructuredGrid  *grid, double time, unsigned int timeStep, bool lastTimeStep)
{
	vtkNew<vtkCPDataDescription> dataDescription;
	dataDescription->AddInput("input");
	dataDescription->SetTimeData(time, timeStep);

	if (lastTimeStep == true)
	{
		// assume that we want to all the pipelines to execute if it
		// is the last time step.
		dataDescription->ForceOutputOn();
	}

	if (Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
	{
		dataDescription->GetInputDescriptionByName("input")->SetGrid(grid);
		Processor->CoProcess(dataDescription.GetPointer());
	}
}


} // end of Catalyst namespace
