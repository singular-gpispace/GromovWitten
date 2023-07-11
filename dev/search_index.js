var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GromovWitten","category":"page"},{"location":"#GromovWitten","page":"Home","title":"GromovWitten","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GromovWitten.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GromovWitten]","category":"page"},{"location":"#GromovWitten.GromovWitten","page":"Home","title":"GromovWitten.GromovWitten","text":"GromovWitten  is a package for computing Gromov-Witten invariant via Feynman Integral  .\n\n\n\n\n\n","category":"module"},{"location":"#GromovWitten.constterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer}","page":"Home","title":"GromovWitten.constterm","text":"constterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, N::Integer)  returns the constant term of the propagator\n\nExamples (without vertex contribution)\n\njulia> constterm(x[1],x[2],3)\n\n3x[1]^6 + 2x[1]^5x[2] + x[1]^4x[2]^2\n\nconstterm(x1::QQMPolyRingElem, x2::QQMPolyRingElem, z1::QQMPolyRingElem, z2::QQMPolyRingElem,aa::Integer, N::Integer)\n\nExamples (without vertex contribution)\n\nhere aa=1 is the order for sfunction series and N=sum_n=1^3g-3 a_i where a=\\[a_1,…,a_{3g-3}\\] is a branche type.\n\njulia> constterm(x[1],x[2],z[1],z[2],1,2)\n\n1//18x[1]^4z[1]^2z[2]^2 + 1//3x[1]^4z[1]^2 + 1//3x[1]^4z[2]^2 + 2x[1]^4 \n\n+ 1//576x[1]^3x[2]z[1]^2z[2]^2 + 1//24x[1]^3x[2]z[1]^2 + 1//24x[1]^3x[2]z[2]^2 + x[1]^3*x[2]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_branchtype-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}}","page":"Home","title":"GromovWitten.feynman_integral_branchtype","text":"feynmanintegralbranchtype( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for a specified branch type a for all ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> feynmanintegralbranchtype(x,q,G,a)\n\n256q[2]^2q[3]*q[6]\n\n\n\nfeynmanintegralbranchtype(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> feynmanintegralbranchtype(x,q,z,G,a,aa=1,g=[1,0,0])\n\n115//3*q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_branchtype_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Vector{Int64}}","page":"Home","title":"GromovWitten.feynman_integral_branchtype_order","text":"feynmanintegralbranchtype_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}  ,Ω::Vector{Int64};l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for a specified branch type a for a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegralbranchtype_order(x,q,G,a,Ω)\n\n128q[2]^2q[3]*q[6]\n\nfeynmanintegralbranchtype_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}, Ω::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegralbranchtype_order(x,q,z,G,a,Ω,aa=1,g=[1,0,0])\n\n115//6*q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Integer}","page":"Home","title":"GromovWitten.feynman_integral_degree","text":"feynmanintegraldegree(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,d::Integer;l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for  over all the partitions of the degree d  for all  ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> feynmanintegraldegree(x,q,G,3)\n\n288q[1]^3 + 32q[1]^2q[2] + 32q[1]^2q[3] + 32q[1]^2q[5] + 32q[1]^2q[6] + 8q[1]q[2]q[5] + 8q[1]q[2]q[6] + 8q[1]q[3]q[5] + 8q[1]q[3]q[6] + 24q[2]^3 + 152q[2]^2q[3] + 8q[2]^2q[5] + 8q[2]^2q[6] + 152q[2]q[3]^2 + 32q[2]q[3]q[5] + 32q[2]q[3]q[6] + 32q[2]q[4]^2 + 8q[2]q[4]q[5] + 8q[2]q[4]q[6] + 8q[2]q[5]^2 + 32q[2]q[5]q[6] + 8q[2]q[6]^2 + 24q[3]^3 + 8q[3]^2q[5] + 8q[3]^2q[6] + 32q[3]q[4]^2 + 8q[3]q[4]q[5] + 8q[3]q[4]q[6] + 8q[3]q[5]^2 + 32q[3]q[5]q[6] + 8q[3]q[6]^2 + 288q[4]^3 + 32q[4]^2q[5] + 32q[4]^2q[6] + 24q[5]^3 + 152q[5]^2q[6] + 152q[5]q[6]^2 + 24q[6]^3\n\nfeynmanintegraldegree( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> feynmanintegraldegree(x,q,z,G,3,aa=1,g=[1,0,0])\n\n115//3q[1]^3 + 1//4q[1]^2q[2] + 1//4q[1]^2q[3] + 1//4q[1]q[2]^2 + 1//2q[1]q[2]q[3] + 1//4q[1]q[3]^2 + 115//3q[2]^3 + 1//4q[2]^2q[3] + 1//4q[2]q[3]^2 + 115//3q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Integer}","page":"Home","title":"GromovWitten.feynman_integral_degree_order","text":"feynmanintegraldegree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for  over all the partitions of the degree d  for a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegraldegree_order(x,q,G,Ω,4)\n\n4q[1]^2q[2]q[3] + 4q[1]q[2]^2q[5] + 4q[1]q[2]^2q[6] + 4q[1]q[3]^2q[5] + 4q[1]q[3]^2q[6] + 176q[2]^4 + 496q[2]^3q[3] + 60q[2]^3q[5] + 60q[2]^3q[6] + 788q[2]^2q[3]^2 + 128q[2]^2q[3]q[5] + 128q[2]^2q[3]q[6] + 4q[2]^2q[4]q[5] + 4q[2]^2q[4]q[6] + 16q[2]^2q[5]^2 + 16q[2]^2q[6]^2 + 496q[2]q[3]^3 + 128q[2]q[3]^2q[5] + 128q[2]q[3]^2q[6] + 4q[2]q[3]q[4]^2 + 48q[2]q[3]q[5]^2 + 4q[2]q[3]q[5]q[6] + 48q[2]q[3]q[6]^2 + 176q[3]^4 + 60q[3]^3q[5] + 60q[3]^3q[6] + 4q[3]^2q[4]q[5] + 4q[3]^2q[4]q[6] + 16q[3]^2q[5]^2 + 16q[3]^2q[6]^2\n\nfeynmanintegraldegree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegraldegree_order(x,q,z,G,o,3,aa=1,g=[1,0,0])\n\n1//24q[1]^2q[2] + 1//24q[1]^2q[3] + 1//24q[1]q[2]^2 + 1//12q[1]q[2]q[3] + 1//24q[1]q[3]^2 + 1//24q[2]^2q[3] + 1//24q[2]q[3]^2 + 115//6q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_sum-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Union{Int64, Vector{Int64}}}","page":"Home","title":"GromovWitten.feynman_integral_degree_sum","text":"feynmanintegraldegree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))\n\ncompute the sum of all Feynman Integrals up to a certain degree d for all  ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> feynmanintegraldegree_sum(x,q,G,3)\n\n288q[1]^3 + 32q[1]^2q[2] + 32q[1]^2q[3] + 32q[1]^2q[5] + 32q[1]^2q[6] + 8q[1]^2 + 8q[1]q[2]q[5] + 8q[1]q[2]q[6] + 8q[1]q[3]q[5] + 8q[1]q[3]q[6] + 24q[2]^3 + 152q[2]^2q[3] + 8q[2]^2q[5] + 8q[2]^2q[6] + 152q[2]q[3]^2 + 32q[2]q[3]q[5] + 32q[2]q[3]q[6] + 8q[2]q[3] + 32q[2]q[4]^2 + 8q[2]q[4]q[5] + 8q[2]q[4]q[6] + 8q[2]q[5]^2 + 32q[2]q[5]q[6] + 8q[2]q[6]^2 + 24q[3]^3 + 8q[3]^2q[5] + 8q[3]^2q[6] + 32q[3]q[4]^2 + 8q[3]q[4]q[5] + 8q[3]q[4]q[6] + 8q[3]q[5]^2 + 32q[3]q[5]q[6] + 8q[3]q[6]^2 + 288q[4]^3 + 32q[4]^2q[5] + 32q[4]^2q[6] + 8q[4]^2 + 24q[5]^3 + 152q[5]^2q[6] + 152q[5]q[6]^2 + 8q[5]q[6] + 24q[6]^3\n\nfeynmanintegraldegree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> feynmanintegraldegree_sum(x,q,z,G,3,aa=1,g=[1,0,0])\n\n115//3q[1]^3 + 1//4q[1]^2q[2] + 1//4q[1]^2q[3] + 19//4q[1]^2 + 1//4q[1]q[2]^2 + 1//2q[1]q[2]q[3] + 1//4q[1]q[2] + 1//4q[1]q[3]^2 + 1//4q[1]q[3] + 1//12q[1] + 115//3q[2]^3 + 1//4q[2]^2q[3] + 19//4q[2]^2 + 1//4q[2]q[3]^2 + 1//4q[2]q[3] + 1//12q[2] + 115//3q[3]^3 + 19//4q[3]^2 + 1//12q[3]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_sum_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Union{Int64, Vector{Int64}}}","page":"Home","title":"GromovWitten.feynman_integral_degree_sum_order","text":"feynmanintegraldegreesumorder(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\ncompute the sum of all Feynman Integrals up to a certain degree d with a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegraldegreesumorder(x,q,G,Ω,3) 12q[2]^3 + 76q[2]^2q[3] + 4q[2]^2q[5] + 4q[2]^2q[6] + 76q[2]q[3]^2 + 16q[2]q[3]q[5] + 16q[2]q[3]q[6] + 4q[2]q[3] + 12q[3]^3 + 4q[3]^2q[5] + 4q[3]^2q[6] –––––––––––––\n\nfeynmanintegraldegreesumorder(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegraldegreesumorder(x,q,z,G,o,3,aa=1,g=[1,0,0])\n\n1//24q[1]^2q[2] + 1//24q[1]^2q[3] + 1//24q[1]q[2]^2 + 1//12q[1]q[2]q[3] + 1//24q[1]q[2] + 1//24q[1]q[3]^2 + 1//24q[1]q[3] + 1//24q[2]^2q[3] + 1//24q[2]q[3]^2 + 1//24q[2]q[3] + 115//6q[3]^3 + 19//8q[3]^2 + 1//24q[3]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.filter_term-Tuple{Union{Int64, Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, Vector{Int64}}","page":"Home","title":"GromovWitten.filter_term","text":"filter_term(p::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, s::Vector{Int64})\n\nreplaces all terms of the polynomial p with zero whenever the variables raised to a power of s1 exceed the specified power s.\n\nExamples\n\njulia> p=feynmanintegraldegree(x,q,G,4)\n\n8q[1]^3q[2] + 8q[1]^3q[3] + 8q[1]^3q[4] + 54q[1]^2q[2]^2 + 18q[1]^2q[2]*q[3] \n\n18q[1]^2q[2]q[4] + 54q[1]^2q[3]^2 + 18q[1]^2q[3]q[4] + 54q[1]^2q[4]^2 \n56q[1]q[2]^3 + 6q[1]q[2]^2q[3] + 6q[1]q[2]^2q[4] + 6q[1]q[2]*q[3]^2 \n12q[1]q[2]q[3]q[4] + 6q[1]q[2]q[4]^2 + 56q[1]q[3]^3 + 6q[1]q[3]^2q[4] \n6q[1]q[3]q[4]^2 + 56q[1]*q[4]^3\n\nwe replace all term in p  with q[1]^a*q[2]^b*q[3]^c > q[1]*q[2]*q[3] by zero,this means all power (abc)(111)\n\njulia> filter_term(p,[q[1],q[2],q[3]],[1,1,1])\n\n12q[1]q[2]q[3]q[4] + 6q[1]q[2]q[4]^2 + 6q[1]q[3]q[4]^2 + 56q[1]q[4]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.flip_signature-Tuple{graph, Vector{Int64}, Vector{Int64}}","page":"Home","title":"GromovWitten.flip_signature","text":"flip_signature(G::graph ,p::Vector{Int64},a::Vector{Int64})\n\nLet   Ω=[x1,...,xn] be  a given Order and a  a branche type,flip_signature  returns -1 if xi<xj and O else.   It will return -2 in case the Graph G has a loop. \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.inv_sfunction-Tuple{Nemo.QQMPolyRingElem, Int64}","page":"Home","title":"GromovWitten.inv_sfunction","text":"inv_sfunction(z::QQMPolyRingElem,aa::Int64) returns the inverse sfunction\n\n frac1S(zaa)=fracz2 Sinh(z2)= \n sum_n =  0^aa left( leftbeginarrayll\n     1  if  n = 1\n     - frac- 2^n (- 2 + 2^n)n B_n  (n  = 1  (- 1\n     + n)mod 2 = 1)\n   endarrayright right) z^n \n\nWhere B_n is Bernoulli number and aa rightarrow infty.\n\nExamples\n\njulia> inv_sfunction(x[1],4)\n\n7//5760x[1]^4 - 1//24x[1]^2 + 1\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.loopterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer, Integer}","page":"Home","title":"GromovWitten.loopterm","text":"loopterm( z::QQMPolyRingElem, q::QQMPolyRingElem, aa::Integer, a::Integer)\n\nreturns loop contribution with nonzero genus gi at a vertex i. \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.partition-Tuple{Integer, Integer}","page":"Home","title":"GromovWitten.partition","text":"partition(k::Integer, n::Integer)    \n\nExamples\n\nThis function returns the number of partitions of n into fixed  k parts. \n\njulia> partition(3,4)\n15-element Vector{Vector{Int64}}:\n [4, 0, 0]\n [3, 1, 0]\n [3, 0, 1]\n [2, 2, 0]\n [2, 1, 1]\n [2, 0, 2]\n [1, 3, 0]\n [1, 2, 1]\n [1, 1, 2]\n [1, 0, 3]\n [0, 4, 0]\n [0, 3, 1]\n [0, 2, 2]\n [0, 1, 3]\n [0, 0, 4]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.proterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer, Integer, Integer}","page":"Home","title":"GromovWitten.proterm","text":"proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer, N::Integer)\n\nreturns the non constant term of the propagator\n\nExamples (without vertex contribution)\n\njulia> proterm(x[1],x[2],q[1],1,2)\n\nx[1]^3x[2]q[1] + x[1]x[2]^3q[1]\n\nExamples (with vertex contribution)\n\nproterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem,z1::QQMPolyRingElem, z2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer,aa::Integer, N::Integer)  julia> proterm(x[1],x[2],z[1],z[2],q[1],1,1,2)\n\n1//576x[1]^3x[2]q[1]z[1]^2z[2]^2 + 1//24x[1]^3x[2]q[1]*z[1]^2 \n\n1//24x[1]^3x[2]q[1]z[2]^2 + x[1]^3x[2]q[1] + 1//576x[1]x[2]^3q[1]z[1]^2*z[2]^2 \n1//24x[1]x[2]^3q[1]z[1]^2 + 1//24x[1]x[2]^3q[1]z[2]^2 + x[1]x[2]^3q[1]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.sfunction-Tuple{Nemo.QQMPolyRingElem, Int64}","page":"Home","title":"GromovWitten.sfunction","text":"sfunction(z::QQMPolyRingElem,k::Int64)\n\nNote:The function sfunction(z,k) takes account vertex contributions. \n\nS(z aa) = sum_n = 0^aa dfrac2^- 1 - n (1 + (- 1)^n)\n(n + 1)  z^n = sum_n = 0^aa dfrac2^- 2 n (2 n + 1) \nz^n aa rightarrow infty\n\nExamples\n\njulia> sfunction(x[1],4)\n\n1//92897280z[1]^8 + 1//322560z[1]^6 + 1//1920z[1]^4 + 1//24z[1]^2 \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.signature_and_multiplicities-Tuple{graph, Vector{Int64}}","page":"Home","title":"GromovWitten.signature_and_multiplicities","text":"signature_and_multiplicities( G::graph, a::Vector{Int64})\n\nreturns flip_signature and their multiplicities.\n\nExamples\n\njulia> G=graph(ve)\n\ngraph([(1, 1), (1, 2), (2, 3), (3, 1)])\n\na=[2,0,0,1]\n\njulia> signature_and_multiplicities(G,a)\n\n4-element Vector{Tuple{Int64, Vector{Int64}}}:\n\n(1, [-2, 0, 0, 1])  (2, [-2, -1, 0, 1])  (2, [-2, 0, -1, 1])  (1, [-2, -1, -1, 1])\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.substitute-Tuple{Vector{Nemo.QQMPolyRingElem}, Union{Int64, Nemo.QQMPolyRingElem}}","page":"Home","title":"GromovWitten.substitute","text":"substitute(q::Vector{QQMPolyRingElem},p::Union{QQMPolyRingElem, Int64})\n\nreplace all the variables by the first variable of p.   With x=\\[x_1,x2,x3 \\] and p(x_1,x2,x3), substitute(x,p) returns p(x1,x_1,x_1)\n\njulia> f=x[1]x[2]+x[1]^3x[2]+5x[1]^6-2x[3]x[2] 5x[1]^6 + x[1]^3x[2] + x[1]x[2] - 2x[2]x[3]\n\njulia> substitute(x,f)\n\n5*x[1]^6 + x[1]^4 - x[1]^2\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.sum_of_coeff-Tuple{Nemo.QQMPolyRingElem}","page":"Home","title":"GromovWitten.sum_of_coeff","text":"sumofcoeff(p::QQMPolyRingElem)\n\ncompute the sum of coefficient of the polynomial p. \n\njulia> f=3x[1]^6 + 2x[1]^5x[2] + x[1]^4x[2]^2\n\njulia> sumofcoeff(f)\n\n6\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten-2","page":"Home","title":"GromovWitten","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package GromovWitten computes generating series for tropical Hurwitz numbers of elliptic curves via mirror symmetry and Feynman integrals, and thus, by a correspondence theorem, Hurwitz numbers in the sense algebraic geometry. Generalizations of the method also allow for the computation of Gromov-Witten invariants for ellptic curves, and are also implemented in the package. GromovWitten is based on the computeralgebra system OSCAR and is provided as a package for the Julia programming language.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We assume that Julia is installed in a recent enough version to run OSCAR. Navigate in a terminal to the folder where you want to install the package and pull the package from Github:","category":"page"},{"location":"","page":"Home","title":"Home","text":"git pull https://github.com/singular-gpispace/GromovWitten.git","category":"page"},{"location":"","page":"Home","title":"Home","text":"I the same folder execute the following command:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia --project","category":"page"},{"location":"","page":"Home","title":"Home","text":"This will activate the environment for our package. In Julia install missing packages:","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.instantiate()","category":"page"},{"location":"","page":"Home","title":"Home","text":"and load our package. On the first run this may take some time.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using GromovWitten  ","category":"page"},{"location":"#Example-of-graph-with-vertex-contribution","page":"Home","title":"Example of graph with vertex contribution","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<img width=\"400\" alt=\"image\" src=\"https://github.com/singular-gpispace/GromovWitten/assets/46294807/01eefb54-3db4-49cf-93a5-4ec2218ba9ce\">","category":"page"},{"location":"","page":"Home","title":"Home","text":"To provide an example on how to use our package, we define a graph G from a list of edges:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> ve = [ (1, 2), (2, 3), (3, 1)]  \njulia> G = graph(ve)","category":"page"},{"location":"","page":"Home","title":"Home","text":"We then define a polynomial ring with all variables required by our implementation:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> R,x,q,z=polynomial_ring(G,\"x\",\"q\",\"z\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here, the indexed variables x correspond to the vertices of the graph, the indexed variables y to the edges of the graph, and the indexed variables z again to the vertices of the graph (the latter to be used in the context of Gromov-Witten invariants with non-trivial Psi-classes).","category":"page"},{"location":"","page":"Home","title":"Home","text":"To compute a Feynman iuntegral, we define a partition  a=003  of degree d=3, a fixed order of vertex o=123 and the genus function g=100. The leak in G is L=000 , aa=1 is the order of the sfunction. We have then","category":"page"},{"location":"","page":"Home","title":"Home","text":" julia> feynman_integral_branchtype_order(x,q,z,G,a,o,aa=1,g=[1,0,0])","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Feynman Integral branch type for all ordering with genus function g  is","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> feynman_integral_branchtype(x,q,z,G,a,aa=1,g=[1,0,0])","category":"page"},{"location":"","page":"Home","title":"Home","text":"also we can compute Feynman Integral of degree 4","category":"page"},{"location":"","page":"Home","title":"Home","text":" julia> feynman_integral_degree(x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0])","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally we substitute all q  variables by q_1","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> substitute(feynman_integral_degree(x,q,z,G,3,aa=1,l=[0,0,0],g=[1,0,0]))","category":"page"},{"location":"#Example-of-graph-without-vertex-contribution-and-loop.","page":"Home","title":"Example of graph without vertex contribution and loop.","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<img width=\"400\" alt=\"image\" src=\"https://github.com/singular-gpispace/GromovWitten/assets/46294807/1b45577b-3c92-464f-81f5-57766dcd189e\">","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> G = graph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3,4)] )\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> R, x, q = polynomial_ring(G, \"x\", \"q\")\n(Multivariate polynomial ring in 10 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3],x[4]], QQMPolyRingElem[q[1], q[2],q[3], q[4], q[5], q[6]])","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> a = [0, 2, 1, 0, 0, 1];","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> o=[1,3,4,2];","category":"page"},{"location":"","page":"Home","title":"Home","text":"\njulia>  feynman_integral_branchtype_order(x,q,G,a,o) \n128*q[2]^2*q[3]*q[6]","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> feynman_integral_branchtype(x, q, G, a)  \n256*q[2]^2*q[3]*q[6]","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> f = feynman_integral_degree(x, q, G, 3)\n288*q[1]^3 + 32*q[1]^2*q[2] + 32*q[1]^2*q[3] + 32*q[1]^2*q[5] + 32*q[1]^2*q[6] + 8*q[1]*q[2]*q[5] + 8*q[1]*q[2]*q[6] + 8*q[1]*q[3]*q[5] + 8*q[1]*q[3]*q[6] + 24*q[2]^3 + 152*q[2]^2*q[3] + 8*q[2]^2*q[5] + 8*q[2]^2*q[6] + 152*q[2]*q[3]^2 + 32*q[2]*q[3]*q[5] + 32*q[2]*q[3]*q[6] + 32*q[2]*q[4]^2 + 8*q[2]*q[4]*q[5] + 8*q[2]*q[4]*q[6] + 8*q[2]*q[5]^2 + 32*q[2]*q[5]*q[6] + 8*q[2]*q[6]^2 + 24*q[3]^3 + 8*q[3]^2*q[5] + 8*q[3]^2*q[6] + 32*q[3]*q[4]^2 + 8*q[3]*q[4]*q[5] + 8*q[3]*q[4]*q[6] + 8*q[3]*q[5]^2 + 32*q[3]*q[5]*q[6] + 8*q[3]*q[6]^2 + 288*q[4]^3 + 32*q[4]^2*q[5] + 32*q[4]^2*q[6] + 24*q[5]^3 + 152*q[5]^2*q[6] + 152*q[5]*q[6]^2 + 24*q[6]^3","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>     substitute(q,feynman_integral_degree_sum(x,q,G,8))\n10246144*q[1]^8 + 3294720*q[1]^7 + 886656*q[1]^6 + 182272*q[1]^5 + 25344*q[1]^4 + 1792*q[1]^3 + 32*q[1]^2","category":"page"},{"location":"#Example-of-graph-with-loop.","page":"Home","title":"Example of graph with loop.","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<img width=\"350\" alt=\"image\" src=\"https://github.com/singular-gpispace/GromovWitten/assets/46294807/c94574b9-eaed-47ca-afb0-a53fdd0b64b3\">","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> G=graph([(1, 1), (1, 2), (2, 3), (3, 1)])\ngraph([(1, 1), (1, 2), (2, 3), (3, 1)])","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> R,x,q=polynomial_ring(G,\"x\",\"q\")\n(Multivariate polynomial ring in 7 variables over QQ, QQMPolyRingElem[x[1], x[2], x[3]], QQMPolyRingElem[q[1], q[2], q[3], q[4]])","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> O=[1,2,3]  ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> a=[ 2,  0, 0, 1]\n 4-element Vector{Int64}:\n 2\n 0\n 0\n 1","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> feynman_integral_branchtype_order(x,q,G,a,O)\n3*q[1]^2*q[4]","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> feynman_integral_branchtype(x,q,G,a)  \n6*q[1]^2*q[4]","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> feynman_integral_degree(x,q,G,3)\n6*q[1]^2*q[2] + 6*q[1]^2*q[3] + 6*q[1]^2*q[4] + 18*q[1]*q[2]^2 + 6*q[1]*q[2]*q[3] + 6*q[1]*q[2]*q[4] + 18*q[1]*q[3]^2 + 6*q[1]*q[3]*q[4] + 18*q[1]*q[4]^2","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> substitute(R,x,q,feynman_integral_sum(x,q,G,8))\n20640*q[1]^8 + 9996*q[1]^7 + 4320*q[1]^6 + 1650*q[1]^5 + 456*q[1]^4 + 90*q[1]^3 + 6*q[1]^2","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral","page":"Function","title":"Feynman Integral","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"CurrentModule = GromovWitten\nDocTestSetup = quote\n  using GromovWitten\nend","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"using GromovWitten","category":"page"},{"location":"Feynman Integral/Feynman/#Graph","page":"Function","title":"Graph","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"A Feynman graph is a (non-metrized) graph Γ without ends with n vertices which are labeled x_1     x_n and with labeled edges q_1     q_r. The graph G is represented as a collection of vertices V and edges E. Each edge is a pair (vw) where both v and w are elements of the set of vertices V.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"<img src=\"https://github.com/singular-gpispace/tropicalfeynman/assets/46294807/c5b4b792-6d2f-418f-b38a-21b3c0187a92\" width=\"300\" />","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> ve=[(1, 1), (1, 2), (2, 3), (3, 1)]\n4-element Vector{Tuple{Int64, Int64}}:\n (1, 1)\n (1, 2)\n (2, 3)\n (3, 1)","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"and","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> G=graph(ve)\ngraph([(1, 1), (1, 2), (2, 3), (3, 1)])","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We then define the Rxq=polynomialring(G) from the graph G.  The polynomial ring has $ 5g-5$ variables, consisting of two sets of variables: x_1x_2x_2g-2 and q_1q_2q_3g-3.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia>   R,x,q,z=polynomial_ring(G,\"x\",\"q\",\"z\")\n(Multivariate polynomial ring in 10 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]])","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral-branche-type","page":"Function","title":"Feynman Integral branche type","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"For a given branch type a, we compute the Specific Feynman Integral of the labeled Graph. a is a list of partition of degree d=3 of Gamma.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> a=[2,0,0,1]\n4-element Vector{Int64}:\n 2\n 0\n 0\n 1","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"o","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"is a fixed order of vertices.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> o=[1,2,3]\n3-element Vector{Int64}:\n 1\n 2\n 3","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We compute the Specific Feynman Integral for the Graph G given a fixed vertex ordering o and the partition of degree a. Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia>  feynman_integral_branchtype_order(x,q,z,G,a,o,aa=0,l=[0,0,0],g=[0,0,0])\n3*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"In the case l=[0,0,0], g=[0,0,0] and  $ aa=0$ we can write simply write.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia>  feynman_integral_branchtype_order(x,q,z,G,a,o)\n3*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We compute the Specific Feynman Integral for the Graph G given a fixed partition of degree a for all vertex ordering o. Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> feynman_integral_branchtype(x,q,z,G,a)\n6*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral-2","page":"Function","title":"Feynman Integral","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We compute the  Feynman Integral of the graph G over all  partitions of the degree d=3  for a fixed ordering o.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> feynman_integral_degree_order(x,q,z,G,o,3) # here d=3\n3*q[1]^2*q[4] + q[1]*q[2]*q[3] + q[1]*q[2]*q[4] + q[1]*q[3]*q[4] + 9*q[1]*q[4]^2\n","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We compute the Feynman integral over all the partitions of the degree d of graph G for all vertex ordering.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia>  feynman_integral_degree(x,q,z,G,3) # here d=3\n6*q[1]^2*q[2] + 6*q[1]^2*q[3] + 6*q[1]^2*q[4] + 18*q[1]*q[2]^2 + 6*q[1]*q[2]*q[3] + 6*q[1]*q[2]*q[4] + 18*q[1]*q[3]^2 + 6*q[1]*q[3]*q[4] + 18*q[1]*q[4]^2","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"We compute the sum of coefficients.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Function","title":"Function","text":"julia> sum_of_coeff( feynman_integral_degree(x,q,z,G,3)\n90","category":"page"}]
}
