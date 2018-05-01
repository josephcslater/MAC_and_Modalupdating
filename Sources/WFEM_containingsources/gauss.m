function [pg,wg]=gauss(n)
% GAUSS provides the locations and weights for Gauss Legendre Integration
%  [p,w]=GAUSS(n) where p are the locations, w the weights, and n is a
%       vector of the number of points in each direction.
%       p is n by length(n), w is n by 1, allowing a single loop to
%       perform the entire integration in your routines.
%
%       Current limitations:
%       Values in n must be <11


% Copyright Joseph C. Slater, 7/29/2002
% Equations from Cook, Malkus, Plesha Third Edition
% Values from Zienkkiewicz 3rd Edition (The Finite Element Method)

% p is the local single dimensional gauss point location
% w is the local single dimensional gauss point weight
% pg is the global gauss point location matrix
% wg is the global gauss point location vector.
ln=length(n);
if ln==2
	n(3)=1;
end

p=[];w=[];
if max(n)>10
	disp('Only n < 11 available. Are you kidding! Aren''t 10 enough!')
else
	m=length(n);
	for i=1:m
		nc=n(i);
		
%% Setting Weights
		if nc==1
			p=0;w=2;
		elseif nc==2
			p=[1;-1]/sqrt(3);w=[1;1];
		elseif nc==3
			p=[sqrt(.6);0;-sqrt(.6)];w=[5;8;5]/9;
		elseif nc==4
			p=[sqrt(3+2*sqrt(1.2));
				sqrt(3-2*sqrt(1.2));
				-sqrt(3-2*sqrt(1.2));
				-sqrt(3+2*sqrt(1.2))]/sqrt(7);
			w=[1;1;1;1]/2-[1;-1;-1;1]/6/sqrt(1.2);
		elseif nc==5
			p=[0.906179845938664;
				0.538469310105683;
				0.0;
				-0.538469310105683;
				-0.906179845938664];
			w=[0.236926885056189;
				0.478628670499366;
				0.568888888888889;
				0.478628670499366;
				0.236926885056189];
		elseif nc==6
			p=[0.932469514203152;
				0.661209386466265;
				0.238619186083197;
				-0.238619186083197;
				-0.661209386466265;
				-0.932469514203152];
			w=[0.171324492379170;
				0.360761573048139;
				0.467913934572691;
				0.467913934572691;
				0.360761573048139;
				0.171324492379170];
		elseif nc==7
			p=[0.949107912342759;
				0.741521185599394;
				0.405845151377397;
				0.0;
				-0.405845151377397;
				-0.741521185599394;
				-0.949107912342759];
			w=[0.129484966168870;
				0.279705391489277;
				0.381830050505119;
				0.417959183673469;
				0.381830050505119;
				0.279705391489277;
				0.129484966168870];
		elseif nc==8
			p=[0.960289856497536;
				0.796666477413627;
				0.525532409916329;
				0.183434642495650;
				-0.183434642495650;
				-0.525532409916329;
				-0.796666477413627;
				-0.960289856497536];
			w=[0.101228536290376;
				0.222381034453374;
				0.313706645877887;
				0.362683783378362;
				0.362683783378362;
				0.313706645877887;
				0.222381034453374;
				0.101228536290376];
		elseif nc==9
			p=[0.968160239507626;
				0.836031107326636;
				0.613371432700590;
				0.324253423403809;
				0.0;
				-0.324253423403809;
				-0.613371432700590;
				-0.836031107326636;
				-0.968160239507626];
			w=[0.081274388361574;
				0.180648160694857;
				0.260610696402935;
				0.312347077040003;
				0.330239355001260;
				0.312347077040003;
				0.260610696402935;
				0.180648160694857;
				0.081274388361574];
		elseif nc==10
			p=[0.973906528517172;
				0.865063366688985;
				0.679409568299024;
				0.433395394129247;
				0.148874338981631;
				-0.148874338981631;
				-0.433395394129247;
				-0.679409568299024;
				-0.865063366688985;
				-0.973906528517172];
			w=[0.066671344308688;
				0.149451349150581;
				0.219086362515982;
				0.269266719309996;
				0.295524224714753;
				0.295524224714753;
				0.269266719309996;
				0.219086362515982;
				0.149451349150581;
				0.066671344308688];
		end
%% Combine information into matrices
   		if (m>1)&&(i==1)
			if m==3
				for j=1:n(1)
					pg((j-1)*n(2)*n(3)+1:j*n(2)*n(3),1)=ones(n(2)*n(3),1)*p(j);
					wg((j-1)*n(2)*n(3)+1:j*n(2)*n(3),1)=ones(n(2)*n(3),1)*w(j);
				end
			else
				%Doesn't work.
				for j=1:n(1)
					pg((j-1)*n(2)+1:j*n(2),1)=ones(n(2),1)*p(j);
					wg((j-1)*n(2)+1:j*n(2),1)=ones(n(2),1)*w(j);
				end
			end

		elseif (m>1)&&(i==2)
			%disp('Now i=2')
			
			if m==3
				for j=1:n(2)
					for k=1:n(3)
						%[(j-1)*n(3)+k:n(2)*n(3):n(1)*n(2)*n(3)]
						pg((j-1)*n(3)+k:n(2)*n(3):n(1)*n(2)*n(3),2)=ones(n(1),1)*p(j);
						wg((j-1)*n(3)+k:n(2)*n(3):n(1)*n(2)*n(3),2)=ones(n(1),1)*w(j);
					end
				end
			else
				for j=1:n(2)
					for k=1:n(3)
						%[(j-1)*n(3)+k:n(2)*n(3):n(1)*n(2)*n(3)]
						pg((j-1)+k:n(2):n(1)*n(2),2)=ones(n(1),1)*p(j);
						wg((j-1)+k:n(2):n(1)*n(2),2)=ones(n(1),1)*w(j);
					end
				end

				wg=wg(:,1).*wg(:,2);
			end

			
			%pg
		elseif (m>1)&&(i==3)
			for j=1:n(3)
				pg(j:n(3):n(1)*n(2)*n(3),3)=ones(n(2)*n(1),1)*p(j);
				wg(j:n(3):n(1)*n(2)*n(3),3)=ones(n(2)*n(1),1)*w(j);
			end
			wg=wg(:,1).*wg(:,2).*wg(:,3);
			%pg
		else
			pg=p;
			wg=w;
		end
	end
end
if ln==2
	wg=wg/2;
	pg=pg(:,1:2);
end

	

%disp('Error. Gauss only implemented for 1-D elements so far.')
