function confidence = confCompute(d_cost, h, w)

confidence = zeros(h, w);
for x = 1:w
    for y = 1:h
        response = squeeze(d_cost(y,x,:));         
        [~,~,response_sorted_peaks] = extrema(response);
        if numel(response_sorted_peaks) <= 1
            response_sorted = sort(response,'ascend');
            confidence(y,x) = response_sorted(2)/response_sorted(1);
        else
            confidence(y,x) = response_sorted_peaks(2)/response_sorted_peaks(1);
        end
    end
end
confidence(~isfinite(confidence)) = 1;
confidence = confidence - 1;