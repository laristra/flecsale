-- set some global variables
e0 = 0.106384
gamma = 1.4

-- mesh sizes
num_cells = {20, 20, 20}
length = {1.2, 1.2, 1.2}

-- compute the reference cell size
delta_vol = 1
delta_r = 1.e-12
for i = 1,3 do
  local dx = length[i]/num_cells[i]
  delta_vol = delta_vol * dx
  delta_r = delta_r + math.pow( dx/2, 2 )
end
delta_r = math.sqrt( delta_r )

-- Begin Main Input
hydro = {

  -- The case prefix and postfixes
  prefix = "sedov_3d",
  postfix = "dat",
  -- The frequency of outputs
  output_freq = "10",
  -- The time stepping parameters
  final_time = 1.0,
  max_steps = 10,
  initial_time_step = 1.e-5,
  CFL = { accoustic = 0.25, volume = 0.1, growth = 1.01 },

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

  -- the equation of state
  eos = {
    type = "ideal_gas",
    gas_constant = 1.4,
    specific_heat = 1.0
  },

  -- the initial conditions
  -- return density, velocity, pressure
  ics = function (x,y,z,t)
    local d = 1.0
    local v = {0,0,0}
    local p = 1.e-6
    local r = math.sqrt( x*x + y*y + z*z )
    if r < delta_r then 
      p = (gamma - 1) * d * e0 / delta_vol
    end
    return d, v, p
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
