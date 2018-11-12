function LagrangeBase(n::Integer)
    if n==1
            [x->1-x, x->x],
            [x->-1, x->1]
    elseif n==2
            [x->2 * (x - 1/2) * (x - 1), x->-4 * x * (x - 1), x->2 * x * (x - 1/2)],
            [x->4 * x - 3, x->-8 * x + 4, x->4 * x - 1]
    elseif n==3
            [Poly([1,-11/2,9,-9/2]),Poly([0,9,-45/2,27/2]),Poly([0,-9/2,18,-27/2]),Poly([0,1,-9/2,9/2])],
            [Poly([-11/2,18,-27/2]),Poly([9,-45,81/2]),Poly([-9/2,36,-81/2]),Poly([1,-9,27/2])]
    elseif n==4
            [Poly([1,-25/3,70/3,-80/3,32/3]),
            Poly([0,16,-208/3,96,-128/3]),
            Poly([0,-12,76,-128,64]),
            Poly([0,16/3,-112/3,224/3,-128/3]),
            Poly([0,-1,22/3,-16,32/3])],
            [Poly([-25/3,140/3,-80,128/3]),
            Poly([16,-416/3,288,-512/3]),
            Poly([-12,152,-384,256]),
            Poly([16/3,-224/3,224,-512/3]),
            Poly([-1,44/3,-48,128/3])]
    elseif n==5
            [Poly([1,-137/12,375/8,-2125/24,625/8,-625/24]),
            Poly([0,25,-1925/12,8875/24,-4375/12,3125/24]),
            Poly([0,-25,2675/12,-7375/12,8125/12,-3125/12]),
            Poly([0,50/3,-325/2,6125/12,-625,3125/12]),
            Poly([0,-25/4,1525/24,-5125/24,6875/24,-3125/24]),
            Poly([0,1,-125/12,875/24,-625/12,625/24])],
            [Poly([-137/12,375/4,-2125/8,625/2,-3125/24]),
            Poly([25,-1925/6,8875/8,-4375/3,15625/24]),
            Poly([-25,2675/6,-7375/4,8125/3,-15625/12]),
            Poly([50/3,-325,6125/4,-2500,15625/12]),
            Poly([-25/4,1525/12,-5125/8,6875/6,-15625/24]),
            Poly([1,-125/6,875/8,-625/3,3125/24])]
    elseif n==6
            [Poly([1,-147/10,406/5,-441/2,315,-1134/5,324/5]),
            Poly([0,36,-1566/5,1044,-1674,1296,-1944/5]),
            Poly([0,-45,1053/2,-4149/2,3699,-3078,972]),
            Poly([0,40,-508,2232,-4356,3888,-1296]),
            Poly([0,-45/2,297,-2763/2,2889,-2754,972]),
            Poly([0,36/5,-486/5,468,-1026,5184/5,-1944/5]),
            Poly([0,-1,137/10,-135/2,153,-162,324/5])],
            [Poly([-147/10,812/5,-1323/2,1260,-1134,1944/5]),
            Poly([36,-3132/5,3132,-6696,6480,-11664/5]),
            Poly([-45,1053,-12447/2,14796,-15390,5832]),
            Poly([40,-1016,6696,-17424,19440,-7776]),
            Poly([-45/2,594,-8289/2,11556,-13770,5832]),
            Poly([36/5,-972/5,1404,-4104,5184,-11664/5]),
            Poly([-1,137/5,-405/2,612,-810,1944/5])]
     elseif n==7
             [Poly([1,-363/20,22981/180,-331681/720,16807/18,-386561/360,117649/180,-117649/720]),
             Poly([0,49,-10927/20,109417/45,-88837/16,991613/144,-352947/80,823543/720]),
             Poly([0,-147/2,43071/40,-1347647/240,170471/12,-151263/8,1529437/120,-823543/240]),
             Poly([0,245/3,-46501/36,133427/18,-2926819/144,4151329/144,-2941225/144,823543/144]),
             Poly([0,-245/4,2009/2,-872935/144,52822/3,-1899191/72,117649/6,-823543/144]),
             Poly([0,147/5,-9849/20,45962/15,-444185/48,1159683/80,-2705927/240,823543/240]),
             Poly([0,-49/6,49931/360,-634207/720,98441/36,-319333/72,1294139/360,-823543/720]),
             Poly([0,1,-343/20,9947/90,-16807/48,84035/144,-117649/240,117649/720])],
             [Poly([-363/20,22981/90,-331681/240,33614/9,-386561/72,117649/30,-823543/720]),
             Poly([49,-10927/10,109417/15,-88837/4,4958065/144,-1058841/40,5764801/720]),
             Poly([-147/2,43071/20,-1347647/80,170471/3,-756315/8,1529437/20,-5764801/240]),
             Poly([245/3,-46501/18,133427/6,-2926819/36,20756645/144,-2941225/24,5764801/144]),
             Poly([-245/4,2009,-872935/48,211288/3,-9495955/72,117649,-5764801/144]),
             Poly([147/5,-9849/10,45962/5,-444185/12,1159683/16,-2705927/40,5764801/240]),
             Poly([-49/6,49931/180,-634207/240,98441/9,-1596665/72,1294139/60,-5764801/720]),
             Poly([1,-343/10,9947/30,-16807/12,420175/144,-117649/40,823543/720])]
     elseif n==8
            [Poly([1,-761/35,59062/315,-4272/5,34208/15,-18432/5,53248/15,-65536/35,131072/315]),
            Poly([0,64,-30784/35,44672/9,-673792/45,235520/9,-1196032/45,131072/9,-1048576/315]),
            Poly([0,-112,9936/5,-587296/45,1956992/45,-733184/9,3915776/45,-2228224/45,524288/45]),
            Poly([0,448/3,-128192/45,102016/5,-1097728/15,145408,-2441216/15,1441792/15,-1048576/45]),
            Poly([0,-140,2764,-186496/9,703552/9,-1466368/9,1712128/9,-1048576/9,262144/9]),
            Poly([0,448/5,-9024/5,626048/45,-2443264/45,5285888/45,-6406144/45,4063232/45,-1048576/45]),
            Poly([0,-112/3,34288/45,-5984,358784/15,-53248,999424/15,-131072/3,524288/45]),
            Poly([0,64/7,-6592/35,67456/45,-274432/45,124928/9,-802816/45,3801088/315,-1048576/315]),
            Poly([0,-1,726/35,-7504/45,30944/45,-14336/9,94208/45,-65536/45,131072/315])],
            [Poly([-761/35,118124/315,-12816/5,136832/15,-18432,106496/5,-65536/5,1048576/315]),
            Poly([64,-61568/35,44672/3,-2695168/45,1177600/9,-2392064/15,917504/9,-8388608/315]),
            Poly([-112,19872/5,-587296/15,7827968/45,-3665920/9,7831552/15,-15597568/45,4194304/45]),
            Poly([448/3,-256384/45,306048/5,-4390912/15,727040,-4882432/5,10092544/15,-8388608/45]),
            Poly([-140,5528,-186496/3,2814208/9,-7331840/9,3424256/3,-7340032/9,2097152/9]),
            Poly([448/5,-18048/5,626048/15,-9773056/45,5285888/9,-12812288/15,28442624/45,-8388608/45]),
            Poly([-112/3,68576/45,-17952,1435136/15,-266240,1998848/5,-917504/3,4194304/45]),
            Poly([64/7,-13184/35,67456/15,-1097728/45,624640/9,-1605632/15,3801088/45,-8388608/315]),
            Poly([-1,1452/35,-7504/15,123776/45,-71680/9,188416/15,-458752/45,1048576/315])]
    end
end

struct Lagrange
    l
    dl
    degree
    function Lagrange(degree::Integer)
        l,dl=LagrangeBase(degree)
        new(l,dl,degree)
    end
end

function Quadrature(E::Lagrange,nq)
    Q=Quad(nq)
    n=E.degree+1
    L=Vector{Vector{Float64}}(undef,n)
    dL=Vector{Vector{Float64}}(undef,n)
    for i = 1:n
        L[i]=E.l[i].(Q.Points)
        dL[i]=E.dl[i].(Q.Points)
    end
    L, dL, Q
end
