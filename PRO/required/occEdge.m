function occEdge = occEdge(occ, im_edge, thre)

output = occ .* im_edge;
[edgelist, edgeim] = edgelink(im_edge);        
for k = 1:size(edgelist,2)
    n = sum(sum(edgeim == k));
    p = sum(sum(occ(edgeim == k)));
    if (p/n > thre)
        output(edgeim == k) = 1;
    end            
end
[~, occEdge] = edgelink(output);
occEdge = occEdge > 0;
