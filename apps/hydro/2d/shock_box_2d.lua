hydro = {
  -- The case prefix and postfixes
  prefix = "shock_box_2d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "7",
  -- The time stepping parameters
  final_time = 0.2,
  max_steps = 1e6,
  CFL = 1./2.,
  -- the mesh
  mesh = {
    type = "box",
    dimensions = {10, 10},
    lengths = {1., 1.}
  },
  -- the initial conditions
  ics = function (x,y)
    if x < 0 and y < 0 then
      return 0.125, {0,0}, 0.1
    else
      return 1.0, {0,0}, 1.0
    end
  end 
}
