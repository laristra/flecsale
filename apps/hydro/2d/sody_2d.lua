hydro = {
  -- The case prefix and postfixes
  prefix = "sody_2d",
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
    dimensions = {1, 100},
    xmin = {-0.05, -0.5},
    xmax = { 0.05,  0.5}
  },
  -- the equation of state
  eos = {
    type = "ideal_gas",
    gas_constant = 1.4,
    specific_heat = 1.0
  },
  -- the initial conditions
  -- return density, velocity, pressure
  ics = function (x,y,t)
    if y < 0 then
      return 1.0, {0,0}, 1.0
    else
      return 0.125, {0,0}, 0.1
    end
  end 
}
