controller= require('ControllerVI')
ffi = require('ffi')
local control= ControllerVI()
--io.write(string.format("'%s' is years old\n", control:learn(ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 1,2)))
c = control:learn(ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)
--c = learn(control,ffi.new("double[3]", {1,2,3}),ffi.new("double[1]", {1}), 3, 1, 1,2)
--print(c[1])

