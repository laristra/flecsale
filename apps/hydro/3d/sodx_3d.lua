hydro = {
  -- The case prefix and postfixes
  prefix = "sodx_3d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "1",
  -- The time stepping parameters
  final_time = 0.2,
  max_steps = 1e6,
  CFL = 1.,
  -- the mesh
  mesh = {
    type = "box",
    dimensions = {100, 1, 1},
    lengths = {1., 0.1, 0.1}
  },
  -- the initial conditions
  -- return density, velocity, pressure
  ics = function (x,y,z,t)
    if x < 0 then
      return 1.0, {0,0,0}, 1.0
    else
      return 0.125, {0,0,0}, 0.1
    end
  end 
}
