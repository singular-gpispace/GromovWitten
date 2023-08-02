ve = [ (1, 2), (2, 3), (3, 1)]
G=graph(ve)
a=[0,0,3];
l=[0,0,0];
Ω=[1,2,3];
o=[1,2,3];
g=[1,0,0];
R, x, q = polynomial_ring(G,"x","q")

@testset "preimg test" begin
    L = [1, 2, 3, 4, 5,3]
    xi = 3

    @test preimg(L, xi) == 3
end
@test partition(0,1)== [[0]]
@test partition(1,0)== [[0]]
p=6*q[1]^2*q[2]^2*q[3]^4 + 12*q[1]^2*q[2]^2*q[3]^2*q[3]^2 + 6*q[1]^2*q[2]^2*q[3]^4 + 56*q[1]^2*q[3]^6 + 6*q[1]^2*q[3]^4*q[3]^2 + 6*q[1]^2*q[3]^2*q[3]^4 + 56*q[1]^2*q[3]^6 + q[1]
@test  filter_term(p,q[1],1)==q[1]
@test  loopterm(z[1],q[1],1,2)==11//192*q[1]^4*z[1]^4 + 3//4*q[1]^4*z[1]^2 + 3*q[1]^4
ve=[]
G=graph(ve)
@test_throws DomainError R,x,q=polynomial_ring(G,"x","q")