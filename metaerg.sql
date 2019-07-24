-- cat metaerg.sql | sqlite3 metaerg.db

CREATE TABLE pfam2go
(
    pfam_acc TEXT,
    go_acc TEXT,
    PRIMARY KEY (pfam_acc, go_acc)
);
CREATE INDEX pfam2go_idx_pfam_acc ON pfam2go (pfam_acc);
.separator "\t"
.import pfam2go.sqldb.tsv pfam2go

CREATE TABLE pfams
(
    acc varchar(15) primary key,
    id varchar(50),
    DE TEXT,
    TP varchar(15),
    ML integer
 );

CREATE INDEX pfams_idx_acc ON pfams (acc);
.separator "\t"
.import Pfam-A.hmm.sql.db.tsv pfams

CREATE TABLE tigrfam2go
(
    tigrfam_acc TEXT,
    go_acc TEXT,
    PRIMARY KEY (tigrfam_acc, go_acc)
);
CREATE INDEX tigrfam2go_idx_tigrfam_acc ON tigrfam2go (tigrfam_acc);
.separator "\t"
.import TIGRFAMS2GO.sqldb.tsv tigrfam2go

CREATE TABLE tigrfams
(
    acc varchar(15) primary key,
    id varchar(50),
    DE TEXT,
    IT varchar(15),
    LENG integer,
    gs varchar(50),
    ec varchar(15),
    mainrole text,
    sub1role text
 );

CREATE INDEX tigrfams_idx_acc ON tigrfams (acc);
.separator "\t"
.import TIGRFAMs.hmm.sqldb.tsv tigrfams


CREATE TABLE Gos(
   goid varchar(25) primary key,
   name text
);
CREATE INDEX Gos_idx_goid on Gos (goid);
.separator "\t"
.import go.sqldb.tsv Gos


CREATE TABLE uniprot_spid2annot(
   spid TEXT PRIMARY KEY,
   KO TEXT,
   KEGG TEXT,
   EC TEXT,
   Desc TEXT,
   GO TEXT,
   OS TEXT,
   OC TEXT
);
CREATE INDEX uniprot_spid2annot_idx_spid ON uniprot_spid2annot (spid);
.separator "\t"
.import sprot.sqldb.tsv  uniprot_spid2annot

CREATE TABLE FOAM_ontology(
   L1 TEXT,
   L2 TEXT,
   L3 TEXT,
   L4 TEXT,
   KO TEXT,
   PRIMARY KEY (L1, L2, L3, L4, KO)
);

.separator "\t"
.import FOAM-onto_rel1.uniq.tsv FOAM_ontology

CREATE TABLE genome2taxon
(
    gid TEXT PRIMARY KEY,
    lineage text 
);
CREATE INDEX genome2taxon_idx_gid ON genome2taxon (gid);
.separator "\t"
.import genomedb.tax.tsv genome2taxon
