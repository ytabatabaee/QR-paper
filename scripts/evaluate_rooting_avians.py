from utils import *

if __name__ == '__main__':
    taxa_num = [30, 25, 20, 15, 10, 5]
    sample_counts = [200, 200, 200, 200, 200, 1000]
    seq_length = ['true', '1500', '1000', '500', '250']
    dataset_path = '/home/syt3/scratch/QR-data/'
    tree_types = ['model', 'astral']
    rooted_tree_types = ['qr_exh.1.2.4.modeltree', 'qr_exh.1.2.4.astraltree',
                         'MV.1.5-RAxML_result.astraltree', 'MP.1.5-RAxML_result.astraltree',
                         'MV.1.5-RAxML_result.modeltree', 'MP.1.5-RAxML_result.modeltree',
                         'rd.1.7.0-RAxML_result.modeltree', 'rd.1.7.0-RAxML_result.astraltree',
                         'rd_exh.1.7.0-RAxML_result.modeltree',
                         'mad.2.2-RAxML_result.modeltree.rooted', 'mad.2.2-RAxML_result.astraltree.rooted',
                         'qr_le.1.2.4.modeltree', 'qr_le.1.2.4.astraltree']
    replicates = range(1, 21)

    df = pd.DataFrame(
        columns=["dataset", "rooted-tree-name", "species-tree", "sample-cnt", "replicate", "cd",
                 "acc", "n", "k", "seq-length"])

    for i in range(len(taxa_num)):
        n = taxa_num[i]
        sample_cnt = sample_counts[i]
        seq_length_adjusted = seq_length
        if n > 10:
            seq_length_adjusted = ['true', '1000']
        rooted_tree_types_adjusted = rooted_tree_types
        if n > 10:
            rooted_tree_types_adjusted = [x for x in rooted_tree_types_adjusted if
                                          x not in ['qr_exh.1.2.4.modeltree', 'qr_exh.1.2.4.astraltree',
                                                    'MP.1.5-RAxML_result.astraltree', 'MV.1.5-RAxML_result.astraltree',
                                                    'rd.1.7.0-RAxML_result.astraltree',
                                                    'mad.2.2-RAxML_result.astraltree.rooted', 'qr_le.1.2.4.astraltree']]
        if n > 5:
            rooted_tree_types_adjusted = [x for x in rooted_tree_types_adjusted if
                                          x not in ['rd_exh.1.7.0-RAxML_result.modeltree']]
        for name in rooted_tree_types_adjusted:
            for l in seq_length_adjusted:
                for r in range(1, 21):
                    cd_avg = 0
                    acc_avg = 0
                    print(n, name, l, r)
                    for s in range(sample_cnt):
                        tns = dendropy.TaxonNamespace()
                        true_species_tree_path = dataset_path + 'species_trees/' + str(n) + '_taxon/' + str(
                            s + 1) + '/model-species.tre'
                        estimated_species_tree_path = dataset_path + 'rooted_trees/' + str(
                            n) + '_taxon/' + str(s + 1) + '/' + l + '/R' + str(r) + '/' + name
                        try:
                            true_species_tree = dendropy.Tree.get(path=true_species_tree_path, schema='newick',
                                                                  rooting="force-rooted", taxon_namespace=tns)
                            estimated_species_tree = dendropy.Tree.get(path=estimated_species_tree_path,
                                                                       schema='newick',
                                                                       rooting="force-rooted", taxon_namespace=tns)
                            _, cd = clade_distance(true_species_tree, estimated_species_tree)
                            cd_avg += cd
                            acc_avg += int(cd == 0)
                            df.loc[len(df.index)] = ["Avian", name, name.split('.')[-1], s, r, cd_avg, acc_avg,
                                                     n, 1000, l]
                        except Exception as e:
                            print(str(e), estimated_species_tree_path)
                    cd_avg /= sample_cnt
                    acc_avg /= sample_cnt
    df.to_csv('qr-avian-sim.csv')
