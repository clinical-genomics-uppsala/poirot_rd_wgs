import os


def compile_output_list(wildcards):
    output_files = []
    types = set([unit.type for unit in units.itertuples()])
    for output in output_json:
        if output == "results/{sample}/upd/{sample}.upd_regions.bed":
            for sample in samples[samples.trio_member == "proband"].index:
                proband_trio_id = samples[samples.index == sample].trioid.iloc[0]
                trio_num = samples[samples.trioid == proband_trio_id].shape[0]
                if trio_num != 3:
                    continue
                else:
                    output_files.append(output.format(sample=sample))
        elif output.startswith("results/somalier_trio/"):
            # Only add somalier_trio outputs if config enabled and trios exist
            if config.get("somalier_trio_extract") and has_trio_samples(samples):
                output_files.append(output)
        else:
            output_files += set(
                [
                    output.format(sample=sample, flowcell=flowcell, lane=lane, barcode=barcode, type=unit_type, panel=str_panel)
                    for sample in get_samples(samples)
                    for unit_type in get_unit_types(units, sample)
                    if unit_type in set(output_json[output]["types"])
                    for flowcell in set([u.flowcell for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for barcode in set([u.barcode for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for lane in set([u.lane for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for str_panel in [panel_list.split(".")[0] for panel_list in config["reference"]["str_panels"]]
                ]
            )

    return list(set(output_files))


def generate_copy_code(workflow, output_json):
    code = ""
    for result, values in output_json.items():
        if values["file"] is not None:
            input_file = values["file"]
            output_file = result
            rule_name = values["name"]
            mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
            partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
            threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
            time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
            copy_container = config.get("_copy", {}).get("container", config["default_container"])
            code += f'@workflow.rule(name="{rule_name}")\n'
            code += f'@workflow.input("{input_file}")\n'
            if rule_name == "_copy_reviewer":  # handle rule that has directory as output
                result_file = "{sample}"
                code += f'@workflow.output(directory("{output_file}"))\n'
            else:
                result_file = os.path.basename(output_file)
                code += f'@workflow.output("{output_file}")\n'
            code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
            code += f'@workflow.container("{copy_container}")\n'
            code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
            code += '@workflow.shellcmd("cp -r {input} {output}")\n\n'
            code += "@workflow.run\n"
            code += (
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
                "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
                "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):\n"
                '\tshell ( "(cp -r {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
            )

    exec(compile(code, "result_to_copy", "exec"), workflow.globals)


generate_copy_code(workflow, output_json)
