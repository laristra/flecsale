hydro = {
  -- The case prefix and postfixes
  prefix = "shock_box_3d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "11",
  -- The time stepping parameters
  final_time = 0.2,
  max_steps = 1e6,
  CFL = 1./3.,
  -- the mesh
  mesh = {
    type = "box",
    dimensions = {10, 10, 10},
    lengths = {1., 1., 1.}
  },
  -- the initial conditions
  ics = function (x,y,z)
    if x < 0 and y < 0 and z < 0 then
      return 0.125, {0,0,0}, 0.1
    else
      return 1.0, {0,0,0}, 1.0
    end
  end 
}
