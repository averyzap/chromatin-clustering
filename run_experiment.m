function out = run_experiment(exp_fn)
    [~,name] = fileparts(exp_fn);
    s = load(exp_fn);
    rng(s.seed);
    out = s.fn(s.args{:});
    ts = char(datetime('now','Format','yyyy-MM-dd-HH-mm-ss'));
    [~, msg] = mkdir('data');
    [~, msg] = mkdir(fullfile('data',name));
    save(fullfile('data',name,'output.mat'), 'out');
end