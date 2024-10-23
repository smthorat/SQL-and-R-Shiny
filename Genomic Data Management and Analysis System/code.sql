CREATE TABLE Gene_New (
    GeneID INT PRIMARY KEY,
    Symbol VARCHAR(50),
    Description TEXT,
    TaxonomicName VARCHAR(100),
    CommonName VARCHAR(50),
    GeneType VARCHAR(50)
);
SELECT * FROM Gene_New LIMIT 20;

CREATE INDEX idx_symbol ON Gene_New(Symbol);
CREATE INDEX idx_gene_type ON Gene_New(GeneType);

# Develop SQL Queries for Analysis

# Retrieve all genes: List all records from the Gene_New table.
SELECT * FROM Gene_New;

# Find genes by symbol: Search for genes with a specific symbol (e.g., “NAT1”).
SELECT * FROM Gene_New WHERE Symbol = 'NAT1';

# Filter by gene type: Retrieve genes that are of a certain type, such as “Protein_Coding.”
SELECT * FROM Gene_New WHERE GeneType = 'Protein_Coding';

# Count genes by gene type: Get a count of genes grouped by their type.
SELECT GeneType, COUNT(*) AS GeneCount
FROM Gene_New
GROUP BY GeneType;

## Advanced Queries

# Find genes containing specific keywords in the description: Search for genes where the description contains a specific term (e.g., “transferase”).
SELECT * FROM Gene_New
WHERE Description LIKE '%transferase%';

# Get gene statistics by taxonomic name: Calculate the number of genes and the distribution of gene types for each species.
SELECT TaxonomicName, GeneType, COUNT(*) AS GeneCount
FROM Gene_New
GROUP BY TaxonomicName, GeneType
ORDER BY TaxonomicName, GeneCount DESC;


# Identify genes with multiple occurrences in the database (if applicable): Check for potential duplicates based on the gene symbol.
SELECT Symbol, COUNT(*) AS Occurrences
FROM Gene_New
GROUP BY Symbol
HAVING COUNT(*) > 1;


# Advanced Database Features
# When we change the gene function it should automatically update the protein table as well.


CREATE TABLE Gene_Log (
    LogID INT AUTO_INCREMENT PRIMARY KEY,
    GeneID INT,
    OldDescription TEXT,
    NewDescription TEXT,
    ChangeDate DATETIME
);


DELIMITER //

CREATE PROCEDURE UpdateGeneDescription(IN p_GeneID INT, IN p_Description TEXT)
BEGIN
    UPDATE Gene_New
    SET Description = p_Description
    WHERE GeneID = p_GeneID;
END//

DELIMITER ;


DELIMITER //

CREATE TRIGGER BeforeGeneUpdate
BEFORE UPDATE ON Gene_New
FOR EACH ROW
BEGIN
    INSERT INTO Gene_Log (GeneID, OldDescription, NewDescription, ChangeDate)
    VALUES (OLD.GeneID, OLD.Description, NEW.Description, NOW());
END//

DELIMITER ;



CREATE VIEW ProteinCodingGenes AS
SELECT GeneID, Symbol, Description, TaxonomicName, CommonName
FROM Gene_New
WHERE GeneType = 'Protein_Coding';

# Lets try the procedure to check if its working or not
# Testing the Stored Procedure
CALL UpdateGeneDescription(1, 'New Description for Gene 1');

# Verify the update
SELECT * FROM Gene_New WHERE GeneID = 1;

# Check the log table
SELECT * FROM Gene_Log WHERE GeneID = 1 ORDER BY ChangeDate DESC;

# Testing the trigger 
UPDATE Gene_New SET Description = 'Another Description' WHERE GeneID = 1;

# Verify Logging
SELECT * FROM Gene_Log WHERE GeneID = 1 ORDER BY ChangeDate DESC;

# Testing the view 
SELECT * FROM ProteinCodingGenes;

# Verify the data consistency 
SELECT GeneID, Symbol, Description, TaxonomicName, CommonName
FROM Gene_New
WHERE GeneType = 'Protein_Coding';
