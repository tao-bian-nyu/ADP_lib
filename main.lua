require('luaADP/ControllerVI')
ffi = require('ffi')
control= ControllerVI()
K = {0.002}
t = 0 -- time
    
function love.load()
    
    Object = require('ADPgame/classic')
    require('ADPgame/window')
    require('ADPgame/bird')
    require('ADPgame/pipe')

    -- window initalization
    W = 750 -- window size
    H = 600
    window = Window(W, H)
    
    -- bird initalization
    birdY = 300
    
    bird = Bird(birdY)
    
    -- pipe initalization (3*pipeDist) = W
    gap = 100 -- gap between upper pipe and lower
    pipeY = love.math.random(0, H - gap) -- pipe starting position
    pipeDist = 250 -- distance between 2 pipes
    pipeX = W + pipeDist/2 -- initial position 
    pipeW = 50 -- pipe with
    
    -- initialize pipe queue, queue length is 4
    pipeImg = love.graphics.newImage("ADPgame/Figures/brick_grey.png")
	  pipeImg:setWrap('repeat','repeat')
    pipeQuad = love.graphics.newQuad(0, 0, pipeW, H, 32,32)
    
    listOfPipes = {}
    for i= 0, 3 do
        table.insert(listOfPipes,
                     Pipe( pipeX + pipeDist*i,
                           love.math.random(0, H - gap),
                            pipeW, gap))
    end
    
    
end

function love.update(dt)
    -- update bird
    bird:AIupdate(dt) -- ADP
    -- bird:update(dt) -- with ode
    -- bird:update2(dt)
    
    
    -- update pipes
    for i, v in ipairs(listOfPipes) do
        v:update(dt)
    end
    
    -- check collision
    bird:checkCollision()
    
end

function love.draw()
    
    -- draw pipes
    for i, v in ipairs(listOfPipes) do
        v:draw(dt)
    end    
    
    -- draw birds
    bird:draw()
    
    drawOptLine()
    
    -- display parameters
    
    -- love.graphics.print(p, 50, 20)
    love.graphics.setColor(102, 217, 239)
    love.graphics.print(K[ind], 50, 20)
    -- love.graphics.print(t, 50, 20)
    love.graphics.setColor(255, 255, 255)
    -- love.graphics.print(opt_Y, 400, 20)
    
end

function love.keypressed(key)
    if key == 'left' then
        love.load()
    end
end

function drawOptLine()
  -- draw optimal line
    
    love.graphics.setColor(102, 217, 239) -- line color
    l_b = (pipeDist - pipeW)/2 -- left bound
    r_b = pipeW + l_b -- right bound
    
    love.graphics.line(listOfPipes[1].X - l_b, listOfPipes[1].Y+0.5*listOfPipes[1].gap, --A
                       listOfPipes[1].X + r_b, listOfPipes[1].Y+0.5*listOfPipes[1].gap,--B
                       listOfPipes[2].X - l_b, listOfPipes[2].Y+0.5*listOfPipes[2].gap, --C
                       listOfPipes[2].X + r_b, listOfPipes[2].Y+0.5*listOfPipes[2].gap,--D
                       listOfPipes[3].X - l_b, listOfPipes[3].Y+0.5*listOfPipes[3].gap, --E
                       listOfPipes[3].X + r_b, listOfPipes[3].Y+0.5*listOfPipes[3].gap)--F
    
    love.graphics.setColor(255, 255, 255)

end
