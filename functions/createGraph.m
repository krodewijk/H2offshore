[rows, cols] = size(g.lat);
numNodes = rows * cols;

curNode = 1

A = sparse(numNodes, numNodes);


fprintf(1,' \nProgress:    ');
for y = 1:rows
    for x = 1:cols
        % curNode
        % sweep from left to right, up and down around 3 sized square around middle 
        for x_in = (x-1):(x+1)
            for y_in = (y-1):(y+1)
                % check if x and y are within bound
                if(x_in > 0 && y_in > 0 && x_in <= cols && y_in <= rows)
                    endNode = g.intrin2node(x_in, y_in);
                    dist = norm([g.X(y_in, x_in), g.Y(y_in, x_in)] - [g.X(y,x), g.Y(y,x)]);
                    A(curNode, endNode) = dist;
                end
            end
        end
        curNode = curNode + 1;
    end
    fprintf(1,'\b\b\b\b%3.0f%%',(y/rows*100));
end

G = graph(A, "lower", "omitselfloops");