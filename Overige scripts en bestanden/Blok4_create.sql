-- Created by Vertabelo (http://vertabelo.com)
-- Last modification date: 2019-05-15 07:06:42.432

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
    defenition text NULL,
    accession text NOT NULL,
    CONSTRAINT protein_pk PRIMARY KEY (name_id)
);

-- Table: protein_attribute
CREATE TABLE protein_attribute (
    protein_id int NOT NULL,
    seq_id int NOT NULL,
    organism_id int NOT NULL,
    name_id int NOT NULL,
    ident_num int NULL,
    pos_num int NULL,
    gap_num int NULL,
    e_value varchar(20) NULL,
    bit_score float NULL,
    ident_perc float NULL,
    query_cov float NULL,
    CONSTRAINT protein_attribute_pk PRIMARY KEY (protein_id)
);

-- Table: sequence
CREATE TABLE sequence (
    seq_id int NOT NULL,
    sequence text NOT NULL,
    header text NOT NULL,
    score text NOT NULL,
    CONSTRAINT sequence_pk PRIMARY KEY (seq_id)
);

-- foreign keys
-- Reference: protein_attribute_organism_attribute (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_organism_attribute FOREIGN KEY protein_attribute_organism_attribute (organism_id)
    REFERENCES organism (organism_id);

-- Reference: protein_attribute_protein (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_protein FOREIGN KEY protein_attribute_protein (name_id)
    REFERENCES protein (name_id);

-- Reference: protein_attribute_sequence (table: protein_attribute)
ALTER TABLE protein_attribute ADD CONSTRAINT protein_attribute_sequence FOREIGN KEY protein_attribute_sequence (seq_id)
    REFERENCES sequence (seq_id);

-- End of file.

