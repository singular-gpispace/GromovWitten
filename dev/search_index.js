var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = GromovWitten","category":"page"},{"location":"#GromovWitten","page":"Home","title":"GromovWitten","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GromovWitten.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [GromovWitten]","category":"page"},{"location":"#GromovWitten.GromovWitten","page":"Home","title":"GromovWitten.GromovWitten","text":"GromovWitten  is a package for computing Gromov-Witten invariant via Feynman Integral  .\n\n\n\n\n\n","category":"module"},{"location":"#GromovWitten.constterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer}","page":"Home","title":"GromovWitten.constterm","text":"constterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, N::Integer)  returns the constant term of the propagator\n\nExamples (without vertex contribution)\n\njulia> constterm(x[1],x[2],3)\n\n3x[1]^6 + 2x[1]^5x[2] + x[1]^4x[2]^2\n\nconstterm(x1::QQMPolyRingElem, x2::QQMPolyRingElem, z1::QQMPolyRingElem, z2::QQMPolyRingElem,aa::Integer, N::Integer)\n\nExamples (without vertex contribution)\n\nhere aa=1 is the order for sfunction series and N=sum_n=1^3g-3 a_i where a=\\[a_1,…,a_{3g-3}\\] is a branche type.\n\njulia> constterm(x[1],x[2],z[1],z[2],1,2)\n\n1//18x[1]^4z[1]^2z[2]^2 + 1//3x[1]^4z[1]^2 + 1//3x[1]^4z[2]^2 + 2x[1]^4 \n\n+ 1//576x[1]^3x[2]z[1]^2z[2]^2 + 1//24x[1]^3x[2]z[1]^2 + 1//24x[1]^3x[2]z[2]^2 + x[1]^3*x[2]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_branchtype-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}}","page":"Home","title":"GromovWitten.feynman_integral_branchtype","text":"feynmanintegralbranchtype( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for a specified branch type a for all ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> feynmanintegralbranchtype(x,q,G,a)\n\n256q[2]^2q[3]*q[6]\n\n\n\nfeynmanintegralbranchtype(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64} ;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> feynmanintegralbranchtype(x,q,z,G,a,aa=1,g=[1,0,0])\n\n115//3*q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_branchtype_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Vector{Int64}}","page":"Home","title":"GromovWitten.feynman_integral_branchtype_order","text":"feynmanintegralbranchtype_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}  ,Ω::Vector{Int64};l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for a specified branch type a for a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegralbranchtype_order(x,q,G,a,Ω)\n\n128q[2]^2q[3]*q[6]\n\nfeynmanintegralbranchtype_order(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph ,a::Vector{Int64}, Ω::Vector{Int64};aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegralbranchtype_order(x,q,z,G,a,Ω,aa=1,g=[1,0,0])\n\n115//6*q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Integer}","page":"Home","title":"GromovWitten.feynman_integral_degree","text":"feynmanintegraldegree(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,d::Integer;l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for  over all the partitions of the degree d  for all  ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> feynmanintegraldegree(x,q,G,3)\n\n288q[1]^3 + 32q[1]^2q[2] + 32q[1]^2q[3] + 32q[1]^2q[5] + 32q[1]^2q[6] + 8q[1]q[2]q[5] + 8q[1]q[2]q[6] + 8q[1]q[3]q[5] + 8q[1]q[3]q[6] + 24q[2]^3 + 152q[2]^2q[3] + 8q[2]^2q[5] + 8q[2]^2q[6] + 152q[2]q[3]^2 + 32q[2]q[3]q[5] + 32q[2]q[3]q[6] + 32q[2]q[4]^2 + 8q[2]q[4]q[5] + 8q[2]q[4]q[6] + 8q[2]q[5]^2 + 32q[2]q[5]q[6] + 8q[2]q[6]^2 + 24q[3]^3 + 8q[3]^2q[5] + 8q[3]^2q[6] + 32q[3]q[4]^2 + 8q[3]q[4]q[5] + 8q[3]q[4]q[6] + 8q[3]q[5]^2 + 32q[3]q[5]q[6] + 8q[3]q[6]^2 + 288q[4]^3 + 32q[4]^2q[5] + 32q[4]^2q[6] + 24q[5]^3 + 152q[5]^2q[6] + 152q[5]q[6]^2 + 24q[6]^3\n\nfeynmanintegraldegree( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> feynmanintegraldegree(x,q,z,G,3,aa=1,g=[1,0,0])\n\n115//3q[1]^3 + 1//4q[1]^2q[2] + 1//4q[1]^2q[3] + 1//4q[1]q[2]^2 + 1//2q[1]q[2]q[3] + 1//4q[1]q[3]^2 + 115//3q[2]^3 + 1//4q[2]^2q[3] + 1//4q[2]q[3]^2 + 115//3q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Integer}","page":"Home","title":"GromovWitten.feynman_integral_degree_order","text":"feynmanintegraldegree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)))\n\ncompute the Feynman Integral for  over all the partitions of the degree d  for a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegraldegree_order(x,q,G,Ω,4)\n\n4q[1]^2q[2]q[3] + 4q[1]q[2]^2q[5] + 4q[1]q[2]^2q[6] + 4q[1]q[3]^2q[5] + 4q[1]q[3]^2q[6] + 176q[2]^4 + 496q[2]^3q[3] + 60q[2]^3q[5] + 60q[2]^3q[6] + 788q[2]^2q[3]^2 + 128q[2]^2q[3]q[5] + 128q[2]^2q[3]q[6] + 4q[2]^2q[4]q[5] + 4q[2]^2q[4]q[6] + 16q[2]^2q[5]^2 + 16q[2]^2q[6]^2 + 496q[2]q[3]^3 + 128q[2]q[3]^2q[5] + 128q[2]q[3]^2q[6] + 4q[2]q[3]q[4]^2 + 48q[2]q[3]q[5]^2 + 4q[2]q[3]q[5]q[6] + 48q[2]q[3]q[6]^2 + 176q[3]^4 + 60q[3]^3q[5] + 60q[3]^3q[6] + 4q[3]^2q[4]q[5] + 4q[3]^2q[4]q[6] + 16q[3]^2q[5]^2 + 16q[3]^2q[6]^2\n\nfeynmanintegraldegree_order( x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64},d::Integer;aa=0,l=zeros(Int,nv(G)),g=zeros(Int,nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegraldegree_order(x,q,z,G,o,3,aa=1,g=[1,0,0])\n\n1//24q[1]^2q[2] + 1//24q[1]^2q[3] + 1//24q[1]q[2]^2 + 1//12q[1]q[2]q[3] + 1//24q[1]q[3]^2 + 1//24q[2]^2q[3] + 1//24q[2]q[3]^2 + 115//6q[3]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_sum-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Union{Int64, Vector{Int64}}}","page":"Home","title":"GromovWitten.feynman_integral_degree_sum","text":"feynmanintegraldegree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; l=zeros(Int, nv(G)))\n\ncompute the sum of all Feynman Integrals up to a certain degree d for all  ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> feynmanintegraldegree_sum(x,q,G,3)\n\n288q[1]^3 + 32q[1]^2q[2] + 32q[1]^2q[3] + 32q[1]^2q[5] + 32q[1]^2q[6] + 8q[1]^2 + 8q[1]q[2]q[5] + 8q[1]q[2]q[6] + 8q[1]q[3]q[5] + 8q[1]q[3]q[6] + 24q[2]^3 + 152q[2]^2q[3] + 8q[2]^2q[5] + 8q[2]^2q[6] + 152q[2]q[3]^2 + 32q[2]q[3]q[5] + 32q[2]q[3]q[6] + 8q[2]q[3] + 32q[2]q[4]^2 + 8q[2]q[4]q[5] + 8q[2]q[4]q[6] + 8q[2]q[5]^2 + 32q[2]q[5]q[6] + 8q[2]q[6]^2 + 24q[3]^3 + 8q[3]^2q[5] + 8q[3]^2q[6] + 32q[3]q[4]^2 + 8q[3]q[4]q[5] + 8q[3]q[4]q[6] + 8q[3]q[5]^2 + 32q[3]q[5]q[6] + 8q[3]q[6]^2 + 288q[4]^3 + 32q[4]^2q[5] + 32q[4]^2q[6] + 8q[4]^2 + 24q[5]^3 + 152q[5]^2q[6] + 152q[5]q[6]^2 + 8q[5]q[6] + 24q[6]^3\n\nfeynmanintegraldegree_sum(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> feynmanintegraldegree_sum(x,q,z,G,3,aa=1,g=[1,0,0])\n\n115//3q[1]^3 + 1//4q[1]^2q[2] + 1//4q[1]^2q[3] + 19//4q[1]^2 + 1//4q[1]q[2]^2 + 1//2q[1]q[2]q[3] + 1//4q[1]q[2] + 1//4q[1]q[3]^2 + 1//4q[1]q[3] + 1//12q[1] + 115//3q[2]^3 + 1//4q[2]^2q[3] + 19//4q[2]^2 + 1//4q[2]q[3]^2 + 1//4q[2]q[3] + 1//12q[2] + 115//3q[3]^3 + 19//4q[3]^2 + 1//12q[3]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.feynman_integral_degree_sum_order-Tuple{Vector{Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, graph, Vector{Int64}, Union{Int64, Vector{Int64}}}","page":"Home","title":"GromovWitten.feynman_integral_degree_sum_order","text":"feynmanintegraldegreesumorder(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\ncompute the sum of all Feynman Integrals up to a certain degree d with a fixed ordering Ω\n\nExamples (without vertex contribution)\n\njulia> G=graph(ve)\n\ngraph([(1, 3), (1, 2), (1, 2), (2, 4), (3, 4), (3, 4)])\n\njulia> a=[0,2,1,0,0,1];\n\njulia> Ω=[1,3,4,2];\n\njulia> feynmanintegraldegreesumorder(x,q,G,Ω,3) 12q[2]^3 + 76q[2]^2q[3] + 4q[2]^2q[5] + 4q[2]^2q[6] + 76q[2]q[3]^2 + 16q[2]q[3]q[5] + 16q[2]q[3]q[6] + 4q[2]q[3] + 12q[3]^3 + 4q[3]^2q[5] + 4q[3]^2q[6] –––––––––––––\n\nfeynmanintegraldegreesumorder(x::Vector{QQMPolyRingElem}, q::Vector{QQMPolyRingElem},z::Vector{QQMPolyRingElem}, G::graph,o::Vector{Int64}, d::Union{Int64, Vector{Int64}}; aa=0, l=zeros(Int, nv(G)), g=zeros(Int, nv(G)))\n\nExamples (with vertex contribution)\n\njulia> G=graph(ve) graph([(1, 2), (2, 3), (1, 3)])\n\njulia> a=[0,0,3];\n\njulia> Ω=[1,2,3];\n\njulia> feynmanintegraldegreesumorder(x,q,z,G,o,3,aa=1,g=[1,0,0])\n\n1//24q[1]^2q[2] + 1//24q[1]^2q[3] + 1//24q[1]q[2]^2 + 1//12q[1]q[2]q[3] + 1//24q[1]q[2] + 1//24q[1]q[3]^2 + 1//24q[1]q[3] + 1//24q[2]^2q[3] + 1//24q[2]q[3]^2 + 1//24q[2]q[3] + 115//6q[3]^3 + 19//8q[3]^2 + 1//24q[3]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.filter_term-Tuple{Union{Int64, Nemo.QQMPolyRingElem}, Vector{Nemo.QQMPolyRingElem}, Vector{Int64}}","page":"Home","title":"GromovWitten.filter_term","text":"filter_term(p::Union{QQMPolyRingElem, Int64}, variables::Vector{QQMPolyRingElem}, s::Vector{Int64})\n\nreplaces all terms of the polynomial p with zero whenever the variables raised to a power of s1 exceed the specified power s.\n\nExamples\n\njulia> p=feynmanintegraldegree(x,q,G,4)\n\n8q[1]^3q[2] + 8q[1]^3q[3] + 8q[1]^3q[4] + 54q[1]^2q[2]^2 + 18q[1]^2q[2]*q[3] \n\n18q[1]^2q[2]q[4] + 54q[1]^2q[3]^2 + 18q[1]^2q[3]q[4] + 54q[1]^2q[4]^2 \n56q[1]q[2]^3 + 6q[1]q[2]^2q[3] + 6q[1]q[2]^2q[4] + 6q[1]q[2]*q[3]^2 \n12q[1]q[2]q[3]q[4] + 6q[1]q[2]q[4]^2 + 56q[1]q[3]^3 + 6q[1]q[3]^2q[4] \n6q[1]q[3]q[4]^2 + 56q[1]*q[4]^3\n\nwe replace all term in p  with q[1]^a*q[2]^b*q[3]^c > q[1]*q[2]*q[3] by zero,this means all power (abc)(111)\n\njulia> filter_term(p,[q[1],q[2],q[3]],[1,1,1])\n\n12q[1]q[2]q[3]q[4] + 6q[1]q[2]q[4]^2 + 6q[1]q[3]q[4]^2 + 56q[1]q[4]^3\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.flip_signature-Tuple{graph, Vector{Int64}, Vector{Int64}}","page":"Home","title":"GromovWitten.flip_signature","text":"flip_signature(G::graph ,p::Vector{Int64},a::Vector{Int64})\n\nLet   Ω=[x1,...,xn] be  a given Order and a  a branche type,flip_signature  returns -1 if xi<xj and O else.   It will return -2 in case the Graph G has a loop. \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.inv_sfunction-Tuple{Nemo.QQMPolyRingElem, Int64}","page":"Home","title":"GromovWitten.inv_sfunction","text":"inv_sfunction(z::QQMPolyRingElem,aa::Int64) returns the inverse sfunction\n\n frac1S(zaa)=fracz2 Sinh(z2)= \n sum_n =  0^aa left( leftbeginarrayll\n     1  if  n = 1\n     - frac- 2^n (- 2 + 2^n)n B_n  (n  = 1  (- 1\n     + n)mod 2 = 1)\n   endarrayright right) z^n \n\nWhere B_n is Bernoulli number and aa rightarrow infty.\n\nExamples\n\njulia> inv_sfunction(x[1],4)\n\n7//5760x[1]^4 - 1//24x[1]^2 + 1\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.loopterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer, Integer}","page":"Home","title":"GromovWitten.loopterm","text":"loopterm( z::QQMPolyRingElem, q::QQMPolyRingElem, aa::Integer, a::Integer)\n\nreturns loop contribution with nonzero genus gi at a vertex i. \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.partition-Tuple{Integer, Integer}","page":"Home","title":"GromovWitten.partition","text":"partition(k::Integer, n::Integer)    \n\nNote:This function returns the number of partitions of n into fixed  k parts. \n\nExamples\n\njulia> partition(3,4)\n\n15-element Vector{Vector{Int64}}:\n\n[4, 0, 0]  [3, 1, 0]  [3, 0, 1]  [2, 2, 0]  [2, 1, 1]  [2, 0, 2]  [1, 3, 0]  [1, 2, 1]  [1, 1, 2]  [1, 0, 3]  [0, 4, 0]  [0, 3, 1]  [0, 2, 2]  [0, 1, 3]  [0, 0, 4]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.proterm-Tuple{Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Nemo.QQMPolyRingElem, Integer, Integer, Integer}","page":"Home","title":"GromovWitten.proterm","text":"proterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer, N::Integer)\n\nreturns the non constant term of the propagator\n\nExamples (without vertex contribution)\n\njulia> proterm(x[1],x[2],q[1],1,2)\n\nx[1]^3x[2]q[1] + x[1]x[2]^3q[1]\n\nExamples (with vertex contribution)\n\nproterm( x1::QQMPolyRingElem, x2::QQMPolyRingElem,z1::QQMPolyRingElem, z2::QQMPolyRingElem, q::QQMPolyRingElem, a::Integer,aa::Integer, N::Integer)  julia> proterm(x[1],x[2],z[1],z[2],q[1],1,1,2)\n\n1//576x[1]^3x[2]q[1]z[1]^2z[2]^2 + 1//24x[1]^3x[2]q[1]*z[1]^2 \n\n1//24x[1]^3x[2]q[1]z[2]^2 + x[1]^3x[2]q[1] + 1//576x[1]x[2]^3q[1]z[1]^2*z[2]^2 \n1//24x[1]x[2]^3q[1]z[1]^2 + 1//24x[1]x[2]^3q[1]z[2]^2 + x[1]x[2]^3q[1]\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.sfunction-Tuple{Nemo.QQMPolyRingElem, Int64}","page":"Home","title":"GromovWitten.sfunction","text":"sfunction(z::QQMPolyRingElem,k::Int64)\n\nNote:The function sfunction(z,k) takes account vertex contributions. \n\nS(z aa) = sum_n = 0^aa dfrac2^- 1 - n (1 + (- 1)^n)\n(n + 1)  z^n = sum_n = 0^aa dfrac2^- 2 n (2 n + 1) \nz^n aa rightarrow infty\n\nExamples\n\njulia> sfunction(x[1],4)\n\n1//92897280z[1]^8 + 1//322560z[1]^6 + 1//1920z[1]^4 + 1//24z[1]^2 \n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.signature_and_multiplicities-Tuple{graph, Vector{Int64}}","page":"Home","title":"GromovWitten.signature_and_multiplicities","text":"signature_and_multiplicities( G::graph, a::Vector{Int64})\n\nreturns flip_signature and their multiplicities.\n\nExamples\n\njulia> G=graph(ve)\n\ngraph([(1, 1), (1, 2), (2, 3), (3, 1)])\n\na=[2,0,0,1]\n\njulia> signature_and_multiplicities(G,a)\n\n4-element Vector{Tuple{Int64, Vector{Int64}}}:\n\n(1, [-2, 0, 0, 1])  (2, [-2, -1, 0, 1])  (2, [-2, 0, -1, 1])  (1, [-2, -1, -1, 1])\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.substitute-Tuple{Vector{Nemo.QQMPolyRingElem}, Union{Int64, Nemo.QQMPolyRingElem}}","page":"Home","title":"GromovWitten.substitute","text":"substitute(q::Vector{QQMPolyRingElem},p::Union{QQMPolyRingElem, Int64})\n\nreplace all the variables by the first variable of p.   With x=\\[x_1,x2,x3 \\] and p(x_1,x2,x3), substitute(x,p) returns p(x1,x_1,x_1)\n\njulia> f=x[1]x[2]+x[1]^3x[2]+5x[1]^6-2x[3]x[2] 5x[1]^6 + x[1]^3x[2] + x[1]x[2] - 2x[2]x[3]\n\njulia> substitute(x,f)\n\n5*x[1]^6 + x[1]^4 - x[1]^2\n\n\n\n\n\n","category":"method"},{"location":"#GromovWitten.sum_of_coeff-Tuple{Nemo.QQMPolyRingElem}","page":"Home","title":"GromovWitten.sum_of_coeff","text":"sumofcoeff(p::QQMPolyRingElem)\n\ncompute the sum of coefficient of the polynomial p. \n\njulia> f=3x[1]^6 + 2x[1]^5x[2] + x[1]^4x[2]^2\n\njulia> sumofcoeff(f)\n\n6\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"To run:\n- pull the package in your local disk\nIn the local folder where the package is located, type in the terminal.\n- julia --project\n-  using GromovWitten\n-  ve=[(1, 2), (1,2),(2, 4), (1, 3) ,(3, 4),(3,4)]\n-   G=graph(ve)\n-   R,x,q=polynomial_ring(G,\"x\",\"q\")\n-  a=[0,2,1,0,0,1]\n-  feynman_integral_branchtype(x,q,G,a)\n-  feynman_integral_degree(x,q,G ,4)\n-  subt(R,x,q,feynman_integral_degree(x,q,G ,4))\n","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral","page":"Feynman Integral","title":"Feynman Integral","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"CurrentModule = GromovWitten\nDocTestSetup = quote\n  using GromovWitten\nend","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"using GromovWitten","category":"page"},{"location":"Feynman Integral/Feynman/#Graph","page":"Feynman Integral","title":"Graph","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"A Feynman graph is a (non-metrized) graph Γ without ends with n vertices which are labeled x_1     x_n and with labeled edges q_1     q_r. The graph G is represented as a collection of vertices V and edges E. Each edge is a pair (vw) where both v and w are elements of the set of vertices V.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"<img src=\"https://github.com/singular-gpispace/tropicalfeynman/assets/46294807/c5b4b792-6d2f-418f-b38a-21b3c0187a92\" width=\"300\" />","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> ve=[(1, 1), (1, 2), (2, 3), (3, 1)]\n4-element Vector{Tuple{Int64, Int64}}:\n (1, 1)\n (1, 2)\n (2, 3)\n (3, 1)","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"and","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> G=graph(ve)\ngraph([(1, 1), (1, 2), (2, 3), (3, 1)])","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We then define the Rxq=polynomialring(G) from the graph G.  The polynomial ring has $ 5g-5$ variables, consisting of two sets of variables: x_1x_2x_2g-2 and q_1q_2q_3g-3.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia>   R,x,q,z=polynomial_ring(G,\"x\",\"q\",\"z\")\n(Multivariate polynomial ring in 10 variables over QQ, Nemo.QQMPolyRingElem[x[1], x[2], x[3]], Nemo.QQMPolyRingElem[q[1], q[2], q[3], q[4]], Nemo.QQMPolyRingElem[z[1], z[2], z[3]])","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral-branche-type","page":"Feynman Integral","title":"Feynman Integral branche type","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"For a given branch type a, we compute the Specific Feynman Integral of the labeled Graph. a is a list of partition of degree d=3 of Gamma.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> a=[2,0,0,1]\n4-element Vector{Int64}:\n 2\n 0\n 0\n 1","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"o","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"is a fixed order of vertices.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> o=[1,2,3]\n3-element Vector{Int64}:\n 1\n 2\n 3","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We compute the Specific Feynman Integral for the Graph G given a fixed vertex ordering o and the partition of degree a. Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia>  feynman_integral_branchtype_order(x,q,z,G,a,o,aa=0,l=[0,0,0],g=[0,0,0])\n3*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"In the case l=[0,0,0], g=[0,0,0] and  $ aa=0$ we can write simply write.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia>  feynman_integral_branchtype_order(x,q,z,G,a,o)\n3*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We compute the Specific Feynman Integral for the Graph G given a fixed partition of degree a for all vertex ordering o. Here we have the defaults values of the leak vector and the genus function  l=[0,0,0], g=[0,0,0] and we set the order of Sfunction $ aa=0$","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> feynman_integral_branchtype(x,q,z,G,a)\n6*q[1]^2*q[4]","category":"page"},{"location":"Feynman Integral/Feynman/#Feynman-Integral-2","page":"Feynman Integral","title":"Feynman Integral","text":"","category":"section"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We compute the  Feynman Integral of the graph G over all  partitions of the degree d=3  for a fixed ordering o.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> feynman_integral_degree_order(x,q,z,G,o,3) # here d=3\n3*q[1]^2*q[4] + q[1]*q[2]*q[3] + q[1]*q[2]*q[4] + q[1]*q[3]*q[4] + 9*q[1]*q[4]^2\n","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We compute the Feynman integral over all the partitions of the degree d of graph G for all vertex ordering.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia>  feynman_integral_degree(x,q,z,G,3) # here d=3\n6*q[1]^2*q[2] + 6*q[1]^2*q[3] + 6*q[1]^2*q[4] + 18*q[1]*q[2]^2 + 6*q[1]*q[2]*q[3] + 6*q[1]*q[2]*q[4] + 18*q[1]*q[3]^2 + 6*q[1]*q[3]*q[4] + 18*q[1]*q[4]^2","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"We compute the sum of coefficients.","category":"page"},{"location":"Feynman Integral/Feynman/","page":"Feynman Integral","title":"Feynman Integral","text":"julia> subt( feynman_integral_degree(x,q,z,G,3)\n90","category":"page"}]
}
