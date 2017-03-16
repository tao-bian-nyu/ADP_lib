controller= require('ControllerVI')
ffi = require('ffi')
control= ControllerVI()
--c = control:learn(ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)


--A = matrix{{-1,2,1},{0.1,-1,2},{0,0.1,1}}
--B = matrix{{0},{0},{1}}
--x = matrix{{2},{0},{0}}
--K = matrix{{2,0,0}}
K = {0.002}

p = 500

--matrix.print(B*2)
dt = 0.001
t = 0

--for i = 1,5 do
for i = 1,5000000 do
	--u = -matrix.mul(K,x)
	if i==1 then

		u = -K[1]*p - K[2]*v - K[3]*a 
	else
	        u = -K[0]*p - K[1]*v - K[2]*a
	end

	u = u + math.sin(t)

	dp = (-0.1*p + v) * dt
	dv = (-0.1*v + a)*dt
	da =  (-0.1 * a + u)*dt
	p = p + dp
	--print(K[0])
	--print(p)
	--print(v)
	--print(a)
	K = control:learn(ffi.new("double[3]", x),ffi.new("double[1]", u), 3, 1, dt,t)
	--print(Kadp[1])
	--print(Kadp[2])
	--print(Kadp[3])
	i = i+1
	t = t + dt
end
--c = learn(control,ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)
--print(c[1])

