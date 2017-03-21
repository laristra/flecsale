#ifndef FLESCALE_ADAPTOR
#define FLESCALE_ADAPTOR

#include <iostream>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>


namespace FlecsaleAdaptor
{
	void Init(int numScripts, char* scripts[]);
	void Finalize();
	void CoProcess(vtkUnstructuredGrid *grid, double time, unsigned int timeStep, bool lastTimeStep);
}

#endif