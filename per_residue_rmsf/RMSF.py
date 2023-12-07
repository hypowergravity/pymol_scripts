from pymol import cmd, stored
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import pandas as pd

matplotlib.use('Agg')
font = {'weight': 'bold', 'family': 'sans-serif', 'sans-serif': ['Helvetica']}

matplotlib.rc('font', **font)


class PymolRMSF:
    def __init__(self, files):
        self.files = files
        self.load_trajectory()
        self.residue = self.select_residue()
        self.remove_frames()
        self.atom_names, self.atom_coord_A = self.atom_name()
        self.mean_coord,self.atom_coord = self.atom_cord()
        self.RMSF,self.RMSF_mean = self.calculte_rmsf()
        self.max = self.colour_by_RMSD()
        self.df,self.residue_df = self.plt_by_atom()
        self.plt_residue("plot")

    def load_trajectory(self):
        for frame in self.files:
            cmd.load(frame + "/output_001.pdb", frame)
            if frame != self.files[0]:
                cmd.align(frame, self.files[0])

    def select_residue(self):
        cmd.select("new", f"{self.files[0]} and alt 'B'")
        my_dict = {'my_list': []}
        cmd.iterate("new", "my_list.append(resi)", space=my_dict)
        residue = "+".join(list(set(list(my_dict.values())[0])))
        return residue

    def remove_frames(self):
        # cmd.remove(f"all and resi {self.residue}")
        [cmd.remove(f"{x} and not alt 'B'")
         for x in self.files[1:]]
        cmd.save("20_togther.pdb", f"{' or '.join(self.files)}")

    def atom_name(self):
        """:returns list of atom names and alt A coordinates"""
        atom_dict = {}
        atom_coord = {}
        for resi in self.residue.split("+"):
            atom_dict[f"d_{resi}"] = []
            cmd.iterate(f"new and resi {resi} and not hydrogen",
                        f"d_{resi}.append(name)", space=atom_dict)
            atom_coord[f"d_{resi}"] = cmd.get_coords(
                f"{self.files[0]} and resi {resi} and not hydrogen and alt 'A'")

        return atom_dict, atom_coord

    def atom_cord(self):
        atom_coord = {}
        for resi in self.residue.split("+"):
            atom_coord[f"d_{resi}"] = []
            for x in self.files:
                atom_coord[f"d_{resi}"].append(cmd.get_coords(
                    f"{x} and resi {resi} and not hydrogen and alt 'B'"))

        mean_coord = {}
        for k in list(atom_coord.keys()):
            mean_coord[k] = np.mean(np.array([atom_coord[k][x].reshape(-1) for x in range(len(atom_coord[k]))]),
                                    axis=0).reshape(atom_coord[k][0].shape)

        return mean_coord,atom_coord

    def calculte_rmsf(self):
        RMSF = {}
        RMSF_mean = {}
        # for resi in list(self.atom_names.keys()):
        #     RMSF[resi] = np.sqrt(np.sum((self.atom_coord_A[resi] -
        #                                  self.mean_coord[resi]) ** 2, axis=1))

        squared_differnce = {}
        squared_mean = {}
        for resi in list(self.atom_names.keys()):
            result = np.zeros(len(self.atom_coord[resi][0]))
            for x in self.atom_coord[resi]:
                # RMSF[resi] = np.sqrt(np.sum((x -
                #                              self.mean_coord[resi]) ** 2, axis=1))
                result += np.sum((x - self.mean_coord[resi]) ** 2, axis=1)

            squared_differnce[resi] = result
            squared_mean[resi] = result.mean()
            RMSF[resi] = np.sqrt(squared_differnce[resi] / len(self.atom_coord[resi]))
            RMSF_mean[resi] =  np.sqrt(squared_mean[resi] / len(self.atom_coord[resi]))


        return RMSF,RMSF_mean

    def colour_by_RMSD(self):
        cmd.create("pdb", f"{self.files[0]}")
        cmd.alter(f'pdb', f'b=0.0')
        text = "alt 'B'"
        for key, value in self.atom_names.items():
            for i, v in enumerate(value):
                cmd.alter(f'pdb and resi {key.split("_")[1]} and name {v} and {text}',
                          f'b={self.RMSF[key][i]}')

        max = list(self.RMSF.values())[0][0]
        for x in list(self.RMSF.values()):
            max = np.max(x) if np.max(x) > max else max
            max = round(max, 2)

        cmd.save(f"color_RMSF_2_B_max_{round(max,2)}.pdb", "pdb")
        cmd.set("valence", "0")
        cmd.hide("all")
        cmd.show_as("cartoon", "pdb")
        cmd.show(
            "sticks", f"pdb and resi {self.residue} and not hydrogen and alt 'B'")
        cmd.cartoon("putty", "pdb")
        cmd.set("cartoon_putty_scale_min", 0, "pdb")
        cmd.set("cartoon_putty_scale_max", max, "pdb")
        cmd.set("cartoon_putty_transform", 0, "pdb")
        cmd.set("cartoon_putty_radius", 0.1, "pdb")
        cmd.spectrum("b", "rainbow", "pdb")
        cmd.ramp_new("count", "pdb", [0, max], "rainbow")
        cmd.recolor()
        cmd.ray("1024", "1024")
        cmd.png(f"protein_coloured_by_RMSF.png")

        return max

    def plt_residue(self,type="bar"):
        # residue_dict_to_plot = {}
        # df_atom_name = self.df
        # df_atom_name = df_atom_name.fillna(0)
        # for c in self.df.columns:
        #     residue_dict_to_plot[c] = np.mean(df_atom_name[df_atom_name[c] != 0][c])

        residue_df = self.residue_df

        fig, ax = plt.subplots(figsize=(15, 8))
        ax.bar(residue_df["Residue_name"], residue_df["RMSF"]) if type == "bar" else ax.plot(residue_df["Residue_name"], residue_df["RMSF"])
        ax.set_ylabel(r'RMSF $\AA$')
        ax.set_xlabel('Residue')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(f"per_residue_RMSF_{type}.png",dpi=300)
        plt.close()


    @staticmethod
    def get_residue_name(resi):
        my_dict = {"my_list": []}
        cmd.iterate(f"pdb and resi {resi}",
                    "my_list.append(resn)", space=my_dict)
        return list(my_dict.values())[0][0]

    def plt_by_atom(self):
        residue_list = {}
        for k, v in self.atom_names.items():
            dict_atom = dict(
                map(lambda i, j: (i, j), self.atom_names[k], self.RMSF[k]))
            myKeys = list(dict_atom.keys())
            myKeys.sort()
            residue_list[k] = {i: dict_atom[i] for i in myKeys}

        residue_list_sorted = [int(x.split("_")[1])
                               for x in list(residue_list.keys())]
        residue_list_sorted.sort()

        residue_dict = {i: residue_list[f"d_{i}"] for i in residue_list_sorted}

        resn_dict = {}
        for x in list(residue_dict.keys()):
            resn_dict[f"{self.get_residue_name(x)}_{x}"] = residue_dict[x]

        resn_residue = {}

        resi = [int(x.split("_")[1]) for x in list(self.RMSF_mean.keys())]
        resi.sort()

        for x in [f"d_{x}" for x in resi]:
            name = x.split("_")[1]
            resn_residue[f"{self.get_residue_name(name)}_{name}"] = self.RMSF_mean[f"d_{name}"]

        df_residue = pd.DataFrame(resn_residue.items(),columns=["Residue_name","RMSF"])
        df_residue.to_csv("residue_wise.csv")


        df_atom_name = pd.DataFrame(resn_dict)
        df_atom_name.to_csv("atom_wise.csv")
        # df_atom_name = df_atom_name.replace(np.nan, 0)
        out_pdf = f'RMSF_by_hevay_atom.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)
        # figs = plt.figure()
        for c in df_atom_name.columns:
            fig, ax = plt.subplots(figsize=(15, 8))
            plt.subplot()
            df_atom_name = df_atom_name.fillna(0)
            plt.plot(df_atom_name.index[df_atom_name[c] != 0],
                     df_atom_name[df_atom_name[c] != 0][c])
            plt.xticks(rotation=90)
            plt.xlabel(f'Residue {c}', weight='bold')
            plt.ylabel('RMSF', weight='bold')
            plt.title(f'Per residue RMSF for heavy atom', weight='bold')
            # plt.grid()
            # plt.legend()
            # plt.rcParams.update({'font.size': 8})
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()

        pdf.close()
        return df_atom_name, df_residue

def run():
    Files = [x for x in os.listdir() if "frame" in x]
    pdb = PymolRMSF(Files)


if __name__ == "__main__":
    run()
