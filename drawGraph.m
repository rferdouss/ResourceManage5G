function [ G] = drawGraph( A, X, Y, area_x, area_y)
%% Given an adjuncy matrix A, with three technologies, it generates a graph and gives graph as output
%   Detailed explanation goes here

N = size(A,1);

AA = zeros(N,N);
for i=1:3
   AA = A(:,:,i) + AA; 
end

c1 = [];
c2 = [];
c3 = [];

% find color of edges
for i = 1 : N
	for j = 1 : N
        if j >= i
            continue;
        end
        if AA(i,j)==1
            c1 = [[i j]; c1];
        elseif AA(i,j)==2
            c2 = [[i j]; c2];
        elseif AA(i,j)==3
            c3 = [[i j]; c3];
        end
    end
end

%draw graph
G = graph(AA);
h = plot(G, 'XData',X,'YData',Y);
axis([0 area_x 0 area_y]);
 ylabel('Meters', 'FontSize',20);
 xlabel('Meters',  'FontSize',20);
 l = legend('mmWaves','Ethernet');  
%color edges according to number of edges between two nodes
if ~isempty(c1)
    s = c1(:,1)';
    t = c1(:,2)';
    highlight(h,s,t,'EdgeColor','b','LineWidth',2)
end
if ~isempty(c2)
    s = c2(:,1)';
    t = c2(:,2)';
    highlight(h,s,t,'EdgeColor','r', 'LineWidth',2, 'LineStyle', '--')
end
if ~isempty(c3)
    s = c3(:,1)';
    t = c3(:,2)';
    highlight(h,s,t,'EdgeColor','g', 'LineWidth',2, 'LineStyle', ':')
end

end

