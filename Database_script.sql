-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2019-05-09 09:02:19.447

-- tables
-- Table: organism
CREATE TABLE organism (
    organism_id int NOT NULL,
    organism_species text NULL,
    organism_genus text NULL,
    organism_family text NULL,
    CONSTRAINT organism_pk PRIMARY KEY (organism_id)
);

-- Table: protein
CREATE TABLE protein (
    name_id int NOT NULL,
    protien_name text NULL,
    description text NULL,
    accession text NOT NULL,
    CONSTRAINT protein_pk PRIMARY KEY (name_id)
);

-- Table: protein_attribute
CREATE TABLE protein_attribute (
    protein_id int NOT NULL,
    sequence_seq_id int NOT NULL,
    organism_attribute_organism_id int NOT NULL,
    protein_name_id int NOT NULL,
    ident_num int NULL,
    pos_num int NULL,
    gap_num int NULL,
    e_value varchar(10) NULL,
    bit_score double(10,8) NULL,
    CONSTRAINT protein_attribute_pk PRIMARY KEY (protein_id)
);

-- Table: sequence
CREATE TABLE sequence (
    seq_id int NOT NULL,
    sequentie text NOT NULL,
    header text NOT NULL,
    score text NOT NULL,
    CONSTRAINT sequence_pk PRIMARY KEY (seq_id)
);

-- foreign keys
-- Reference: protein_attribute_organism_attribute (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_organism_attribute FOREIGN KEY protein_attribute_organism_attribute (organism_attribute_organism_id)
    REFERENCES organism (organism_id);

-- Reference: protein_attribute_protein (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_protein FOREIGN KEY protein_attribute_protein (protein_name_id)
    REFERENCES protein (name_id);

-- Reference: protein_attribute_sequence (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_sequence FOREIGN KEY protein_attribute_sequence (sequence_seq_id)
    REFERENCES sequence (seq_id);

-- End of file.

