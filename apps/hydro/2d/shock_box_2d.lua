hydro = {
  -- The case prefix and postfixes
  prefix = "shock_box_2d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "1e6",
  -- The time stepping parameters
  final_time = 0.2,
  max_steps = 20,
  CFL = 1./2.,
  -- the equation of state
  eos = {
    type = "ideal_gas",
    gas_constant = 1.4,
    specific_heat = 1.0
  },
  -- the initial conditions
  -- return density, velocity, pressure
  ics = function (x,y,t)
    if x < 0 and y < 0 then
      return 0.125, {0,0}, 0.1
    else
      return 1.0, {0,0}, 1.0
    end
  end 
}
