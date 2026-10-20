function [loss, failed, relative_error] = fit_error(estimate, truth, solver_exit, normalized_residual, cfg)
% Per-scenario loss, averaged over tissues and noise by evaluate_protocol.
relative_error = (estimate(:,1:3)-truth(:,1:3))./truth(:,1:3);
failed = solver_exit(:) <= 0 | ~isfinite(solver_exit(:)) | ...
    any(~isfinite(estimate),2) | ~isfinite(normalized_residual(:)) | ...
    estimate(:,2) < 5 | normalized_residual(:) > cfg.residual_threshold | ...
    any(estimate < cfg.fitopts.lb' | estimate > cfg.fitopts.ub',2);
loss = sum(relative_error.^2,2);
loss(~isfinite(loss)) = cfg.nonfinite_error;
loss = loss + cfg.failure_penalty*double(failed);
end
