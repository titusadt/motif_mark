#!/usr/bin/env python

import cairo
import argparse
import re

def get_args():
    parser = argparse.ArgumentParser(description="A proram to find motifs in a sequence")
    parser.add_argument("-f", "--fasta_file", help="to specify fasta input file", required=True)
    parser.add_argument("-m", "--motif_file", help="to specify the motif input file", required=True)
    parser.add_argument("-o", "--out_file", help="to specify the output image file")#required=True)
    return parser.parse_args()

args=get_args() 
fasta_file:str =args.fasta_file
motif_file:str =args.motif_file
png_file:str =args.out_file


def oneline_fasta(fasta_file):
    '''This function takes in a fast file and returns a dictionary with headers as keys and sequences as values'''
    fasta_dict = {}
    header = ''
    sequence = ''
    with open(fasta_file, 'r') as fa_file:
        for line in fa_file:
            line = line.strip()
            if line.startswith('>'):
                if header != '':
                    fasta_dict[header] = sequence
                    sequence = ''
                header = line
            else:
                sequence += line
        # Add the last entry to the dictionary
        if header != '':
            fasta_dict[header] = sequence
    return fasta_dict
                    
#create dictionary to find ambiguity
iupac = {'A':'A', 'C':'C', 'G':'G', 'T':'T', 'U':'[TU]',
         'W':'[ATU]','S':'[CG]', 'M':'[AC]', 'K':'[GTU]',
         'R':'[AG]', 'Y':'[CTU]','B':'[CGTU]','D':'[AGTU]',
          'H':'[ACTU]','V':'[ACG]', 'N':'[ACGTU]' }


class Gene:
    def __init__(self):
        '''creating the gene'''
        self.gene_name = None

    def get_name(self, key):
        key_vals = key.split()
        self.gene_name = key_vals[0][1:]
        return self.gene_name

class Intron:
    def __init__(self):
        self.color = 0, 0, 0, 1 
        self.intron_length = None
    
    def get_intron_length(self,sequence):
        self.intron_length = len(sequence)
        return self.intron_length
    
class Exon:
    def __init__(self):
        self.color =  0, 0, 0
        self.exon_pos = None
    
    def get_exon_pos(self, sequence):
        self.exon_pos = []
        for base in range(len(sequence)):
            if sequence[base].isupper():
                self.exon_pos.append(base)
        return self.exon_pos

class MotifClass:
    def __init__(self):
        '''motif definitaions'''
        self.color_pallette = [(255/255,0,0),(255/255,51/255,255/255),(128/255,1,0),
                        (0,1,1),(155/255,153/255,51/255),(128/255,128/255,128/255),(1,1,0),(128/255,128/255,0),(1,0,1),(1,20/255,147/255)]
        self.motifs = []
        self.match_positions = []
        self.regex_str = ''
        self.regex_list = []
        self.motif_dictionary = {}
    
    def get_motifs(self, motif_file):
        with open(motif_file, 'r') as motifs_file:
            for line in motifs_file:
                line = line.strip().upper()
                self.motifs.append(line)
        return self.motifs
    
    def get_regex(self):
        for motif in self.motifs:
            motif = motif.upper()
            replaced_motif = []

            for base in motif:
                if base in iupac:
                    replaced_motif.append(iupac[base])
                else:
                    replaced_motif.append(base)

            self.regex_str = ''.join(replaced_motif)
            self.regex_list.append(self.regex_str)
        return self.regex_list
    
    def find_motifs(self, sequence,ctx,vertical_pos):
        sequence = sequence.upper()
        for pattern in range(len(self.regex_list)):
            for match in re.finditer("(?=" + self.regex_list[pattern] + ")", sequence):
                end_position = len(self.motifs[pattern])+match.end()
                for position in range(match.start(), end_position):
                    ctx.set_source_rgb(*self.color_pallette[pattern])
                    ctx.rectangle(position+30, vertical_pos,2,17)
                    ctx.fill()
                # eg key = YGCY, Value  = (start_position(10), end_position(14), color_coord(0,0,1))
                self.motif_dictionary[self.motifs[pattern]] = (match.start(), end_position, self.color_pallette[pattern])
        return self.motif_dictionary
    
                
class Draw:
    def __init__(self):
        '''Picture things'''
        self.WIDTH, self.HEIGHT  = 1000, 100
        self.intron_length = None
        self.text_pos = 20
        self.legend_vert_pos = 300
        self.legend_text_pos = 315
        self.exon_vert_pos = 35
        

    def draw_intron(self, Intron, ctx, sequence, y_coord):
        intron_len = Intron.get_intron_length(sequence)
        ctx.set_source_rgba(*Intron.color) 
        ctx.set_line_width(2)
        ctx.move_to(30, y_coord)
        ctx.line_to((intron_len+30), y_coord)
        ctx.stroke()
    
    def draw_exon(self, Exon, sequence,ctx):
        #getting the exon positions
        exon_positions = Exon.get_exon_pos(sequence)
        for pos in exon_positions:
            ctx.set_source_rgb(*Exon.color)
            #ctx.rectangle(left_pos, vertical_pos, width,length)
            ctx.rectangle(pos+30,self.exon_vert_pos,2,17)        #(x0,y0,x1,y1)
            ctx.fill() 
            
    
    def draw_motif(self, MotifClass,ctx,vertical_pos):
        for key, value in MotifClass.motif_dictionary.items():
            for position in range(value[0], value[1]):
                ctx.set_source_rgb(*value[2])
                ctx.rectangle(position+30, vertical_pos,2,30)
                ctx.fill()


    def draw_text(self, gene_header,ctx):
        ctx.set_source_rgb(0, 0, 0)
        ctx.select_font_face("Arial",
                            cairo.FONT_SLANT_NORMAL,
                            cairo.FONT_WEIGHT_NORMAL)
        ctx.move_to(40, self.text_pos)
        ctx.show_text(gene_header)
    
    def draw_legend(self, MotifClass, ctx):
        for motif in range(len(MotifClass.motifs)):
            ctx.set_source_rgb(*MotifClass.color_pallette[motif])
            ctx.rectangle(900, self.legend_vert_pos, 20, 20)
            ctx.fill()

            ctx.set_source_rgb(0, 0, 0)
            ctx.move_to(920, self.legend_text_pos)
            ctx.show_text(MotifClass.motifs[motif])
            self.legend_vert_pos +=20
            self.legend_text_pos +=20

            
                 
    
def main():
    
    #create class instance    
    gene = Gene()
    intron = Intron()
    draw = Draw()
    exon = Exon()
    motifs = MotifClass()
    
    #getting fasta header and sequence
    fasta = oneline_fasta(fasta_file)

    motifs.get_motifs(motif_file)
    motifs.get_regex()

    WIDTH, height  = 1000, 200
    y_coord =50
    vertical_pos = 50

    surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, WIDTH, height*len(fasta))
    ctx = cairo.Context(surface)
    ctx.set_source_rgb(1, 1, 1)  # White background
    ctx.paint()

    draw.draw_legend(motifs,ctx)

    for key, value in fasta.items():
        motifs.find_motifs(value,ctx,vertical_pos)
        gene_header = gene.get_name(key)

        #drawing the things
        draw.draw_text(gene_header,ctx)
        draw.draw_intron(intron, ctx, value, y_coord)
        draw.draw_exon(exon, value, ctx)

        #changing the position of things
        y_coord += 200
        vertical_pos += 200
        draw.text_pos +=200
        draw.exon_vert_pos +=200

    surface.write_to_png(png_file)

if __name__=="__main__":
    main()