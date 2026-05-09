# placeholder for targeted extraction of features


from asari.tools.file_io import read_features_from_asari_table





def export_targets(parameters, registry, sample_id, outfile):
    # extract targeted m/z features
    if 'target' in self.parameters and self.parameters['target']:  
        matched_list, _, target_unmapped = all_mass_paired_mapping(
            filtered_FeatureTable['mz'].to_list(), self.parameters['target'], self.parameters['mz_tolerance_ppm']
        )
        print("\nIn targeted extraction, %d target mz values are not found in this dataset: " %len(target_unmapped))
        print('    ', [self.parameters['target'][ii] for ii in target_unmapped])
        matched_targets = [self.parameters['target'][ii[1]] for ii in matched_list]
        targeted_table = filtered_FeatureTable.iloc[[x[0] for x in matched_list], :]
        targeted_table.insert(0, "query_target", matched_targets)
        outfile = os.path.join(self.parameters['outdir'], 'targeted_extraction__'+self.parameters['output_feature_table'])
        targeted_table.to_csv(outfile, index=False, sep="\t")
        print("Targeted extraction Feature table (%d x %d) was written to %s.\n" %(
                            targeted_table.shape[0], number_of_samples, outfile))