-- mesh sizes
num_cells = {100, 10, 10}
length = {1., .1, .1}

-- Begin Main Input
hydro = {

  -- The case prefix and postfixes
  prefix = "sodx_3d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "1",
  -- The time stepping parameters
  final_time = 0.2,
  max_steps = 1e6,
  initial_time_step = 1.,
  CFL = { accoustic = 1., volume = 1., growth = 1.e3 },

  -- the mesh
  mesh = {
    type = "box",
    dimensions = num_cells,
    xmin = {0, 0, 0},
    xmax = length
  },

  -- the equation of state
  eos = {
    type = "ideal_gas",
    gas_constant = 1.4,
    specific_heat = 1.0
  },

  -- the initial conditions
  -- return density, velocity, pressure
  ics = function (x,y,z,t)
    if x < 0.5*length[1] then
      return 1.0, {0,0,0}, 1.0
    else
      return 0.125, {0,0,0}, 0.1
    end
  end,

  -- the boundary conditions
  --
  -- - both +ve and -ve side boundaries can be installed at once since 
  --   they will never overlap
  -- - if they did overlap, then you need to define them seperately or else
  --   its hard to count the number of different conditions on points or edges.
  bcs = {

    -- the +/- x-axis boundaries
    [1] = { 
      type = "symmetry", 
      func = function (x,y,z,t)
        if x == 0 or x == length[1] then
          return true
        else
          return false
        end
      end
    },

    -- the +/- y-axis boundaries
    [2] = { 
      type = "symmetry", 
      func = function (x,y,z,t)
        if y == 0 or y == length[2] then
          return true
        else
          return false
        end
      end
    },

    -- the +/- z-axis boundaries
    [3] = { 
      type = "symmetry", 
      func = function (x,y,z,t)
        if z == 0 or z == length[3] then
          return true
        else
          return false
        end
      end
    }

  } -- bcs

} -- hydro
