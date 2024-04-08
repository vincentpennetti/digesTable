from bandwagon import compute_digestion_bands
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio import SeqIO
import itertools
import math

from pydna.fakeseq import FakeSeq
from pydna.gel import gel, interpolator
from pydna.ladders import NEB_1kb_extend, UPbio_1kb
from PIL import Image, ImageDraw, ImageFont
import time
import os

def main():
    combinations = combo_gen(restriction_enzymes, max_enzymes)
    optimal_combination = ranker(combinations, plasmid_seq)

    # Output optimal combination of restriction enzymes
    if optimal_combination:
        print("Optimal combination of restriction enzymes: ", optimal_combination)

        # if you get a NaN integer operation error whatever, remember it is because the custom ladders need to be put
        # into ladders.py and then imported to gel correctly
        plot_bands(optimal_combination, ladder).show()
    else:
        print("No optimal combination of restriction enzymes found.")


# def calculate_relative_position(band_lengths):
#     """
#     Calculates the relative position of bands in a DNA gel electrophoresis based on a logarithmic transformation of the length of the bands.
#
#     Parameters:
#         band_lengths (list of floats): A list of the lengths of the bands in the DNA gel electrophoresis.
#
#     Returns:
#         list of floats: A list of the relative positions of the bands in the DNA gel electrophoresis.
#     """
#     # Apply a logarithmic transformation to the band lengths
#     log_lengths = [math.log(length) for length in band_lengths]
#
#     # Calculate the total length of all the logarithmic band lengths
#     total_log_length = sum(log_lengths)
#
#     # Calculate the relative position of each band as the cumulative sum of the logarithmic lengths of all the previous bands
#     relative_positions = [sum(log_lengths[:i]) / total_log_length for i in range(1, len(log_lengths) + 1)]
#
#     return relative_positions


#####
def eval_band_pos(bands, interpolator):
    gel_length = 600
    margin = 50
    start_gel = 10

    scale = (gel_length - margin) / interpolator(min(interpolator.x))

    band_pos = []
    for band in bands:
        peak_centre = interpolator(band) * scale + start_gel
        band_pos.append(peak_centre)

    for i in range(len(band_pos)-1):
        if band_pos[i] - band_pos[i+1] < min_distance:
            return False

    return len(band_pos)

#####


# Define function to check if a combination of enzymes meets the criteria
def check_combination(combination, plasmid_seq):
    bands = compute_digestion_bands(plasmid_seq, combination, linear=False)

    if min(bands) < min_fragment_size or max(bands) > max_fragment_size or len(bands) < min_fragments or len(bands) > max_fragments:
        return False

    # rel_pos = calculate_relative_position(bands)
    #
    # discrete = [pos*gel_len for pos in rel_pos]
    # discrete.sort()
    #
    # for position in range(len(discrete) - 1):
    #     p1 = discrete[position]
    #     p2 = discrete[position+1]
    #
    #     if p2 - p1 < min_distance:
    #         return False



    #return len(bands)
    return eval_band_pos(bands, interpolator=interpolator(mwstd=ladder))


# a function for ranking the output RE combinations
# need to redefine rank, at present it is meaningless and unused. The sorting calls do the main job
def ranker(combinations, plasmid_seq):
    # Find optimal combination of restriction enzymes
    baseline_combination = []
    for combination in combinations:
        bands = check_combination(combination, plasmid_seq)
        if bands:
            no_bands = bands
            rank = len(plasmid_seq)/1000 - no_bands

            if len(combination) > 2:
                rank = rank - len(combination)

            baseline_combination.append((combination, rank, no_bands))
        else:
            pass

    baseline_combination.sort(key=lambda tup: tup[2], reverse=True)  # sorts based on number of fragments produced
    baseline_combination.sort(key=lambda tup: len(tup[0]), reverse=False)  # sorts based on number of enzymes required

    return baseline_combination


# Generate all possible combinations of restriction enzymes
def combo_gen(restriction_enzymes, n):
    combinations = []
    for i in range(1, n + 1):
        for j in itertools.combinations(restriction_enzymes, i):
            combinations.append(list(j))
    return combinations


def RE_firstpass(parent_enzymes):
    trimmed = []
    for enzyme in parent_enzymes:
        tmp = enzyme.search(Seq(plasmid_seq), linear=False)
        numcuts = len(tmp)
        if numcuts != 0 and numcuts < max_fragments-1:
            trimmed.append(enzyme)
    return trimmed

# function for plotting the bands
def plot_bands(optimal_combinations, ladder):

    ## don't use deseq for computing the bands--it messes up with some enzymes (like draI using ehamet template)
    ## instead, use bandwagon and generate fakeseqs for the gel image production

    #bands = compute_digestion_bands(plasmid_seq, combination, linear=False)
    # deseq = Dseqrecord(plasmid_seq, circular=True)
    lanes = [ladder]
    #
    # for i in range(5):
    #     lanes.append(deseq.cut(optimal_combinations[i][0]))

    for i, item in enumerate(optimal_combinations):
        if i < 5:
            bands = compute_digestion_bands(plasmid_seq, optimal_combinations[i][0], linear=False)
            fake_bands = [FakeSeq(band) for band in bands]
            lanes.append(fake_bands)
        else:
            break

    # for i in range(5):
    #     bands = compute_digestion_bands(plasmid_seq, optimal_combinations[i][0], linear=False)
    #     fake_bands = [FakeSeq(band) for band in bands]
    #     lanes.append(fake_bands)

    gel_image = gel(samples=lanes, interpolator=interpolator(mwstd=ladder))

    #font = ImageFont.truetype("/Users/vpennetti/PycharmProjects/digestable/venv/lib/DejaVuSans.ttf", 16)

    cwd = os.getcwd()
    filename = "/venv/DejaVuSans.ttf"
    full_path = f"{cwd}{filename}"
    print(full_path)
    font = ImageFont.truetype(full_path, 16)

    draw = ImageDraw.Draw(gel_image)

    # width = int(60 + (lane_width + lanesep) * len(samples))
    # lane width is defaulted to 60, lanesep 10, number of samples is 5 for top 5 digests
    gel_width = int(70*len(lanes)+1)

    # defaulted to 600
    gel_length = 600

    annotated_gel = Image.new('RGB', (gel_width+200, gel_length + 50), (0, 0, 0, 0))
    annotated_gel.paste(gel_image, (0,50))

    annotated_draw = ImageDraw.Draw(annotated_gel)

    for i, item in enumerate(optimal_combinations):
        if i < 5:
            annotated_draw.text((130 + i * 60, 25), str(i + 1), (255, 255, 255), font=font)
            annotated_draw.text((60 * (len(lanes)+1), 50 + (25 * i)),
                                "{} -- {} -- {}f".format(i + 1, str(optimal_combinations[i][0])[1:-1],
                                                         str(optimal_combinations[i][2])), (255, 255, 255), font=font)
        else:
            break

    return annotated_gel




if __name__ == '__main__':
    start = time.time()
    # Define input parameters
    fasta_file = "/Users/vpennetti/Downloads/plasmid1.fasta"

    # Read in plasmid sequence from FASTA file
    with open(fasta_file, "r") as f:
        plasmid_record = SeqIO.read(f, "fasta")
    plasmid_seq = str(plasmid_record.seq)

    parrott_RE = [AccI, AccIII, AclI,AflII,AflIII,AgeI,ApaI,ApaLI,AscI,AvaI,AvrII,BamHI,BglII,BlpI,BspHI,BsrGI,
                  BstYI,BstZ17I,ClaI,DpnI,DraI,EcoO109I,EcoRI,EcoRV,HaeII,HindII,HindIII,HpaI,KasI,Kpn2I,KpnI,MluI,
                  MscI,MspI,NaeI,NcoI,NdeI,NgoMIV,NheI,NotI,PacI,PmeI,PmlI,PsiI,PstI,PvuII,RsaI,SacI,SacII,SalI,SbfI,
                  ScaI,SmaI,SnaBI,SpeI,SphI, SspI, StuI, StyI, SwaI, XbaI, XhoI, XmnI, AarI, BsaI, Esp3I]

    subset_RE = [BamHI, BsrGI, EcoRI, EcoRV, HindIII, KpnI, MluI, Esp3I, NcoI, NdeI, NheI, NotI, PacI, PmeI, SacI,
                 SalI, SbfI, SmaI, XbaI, XhoI]

    min_fragment_size = 500  # bp
    max_fragment_size = 18000  # bp
    min_distance = 7  # mm
    min_fragments = 3
    max_fragments = 9
    max_enzymes = 3
    gel_len = 40

    ## NEB_1kb_extend or UPbio_1kb--remember that they had to be added to the ladders.py file

    if len(plasmid_seq) > 10500:
        ladder = NEB_1kb_extend
    else:
        ladder = UPbio_1kb

    restriction_enzymes = RE_firstpass(parrott_RE)

    main()

    end =time.time()

    print(end-start)
