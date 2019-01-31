package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2013 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.intermine.metadata.Model;
import org.intermine.model.FastPathObject;
import org.intermine.model.InterMineObject;
import org.intermine.model.bio.BioEntity;
import org.intermine.model.bio.DataSet;
import org.intermine.model.bio.Location;
import org.intermine.model.bio.Organism;
import org.intermine.model.bio.SequenceFeature;
import org.intermine.objectstore.ObjectStore;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.util.DynamicUtil;

/**
 * A fasta loader that understand the headers of Medicago fasta CDS fasta files and can make the
 * appropriate extra objects and references.
 * @author Kim Rutherford
 * @author Vivek Krishnakumar
 */
public class MedicagoCDSFastaLoaderTask extends MedicagoFeatureFastaLoaderTask
{
   /**
     * {@inheritDoc}
     */
    @Override
    protected void extraProcessing(Sequence bioJavaSequence,
            org.intermine.model.bio.Sequence flymineSequence,
            BioEntity bioEntity, Organism organism, DataSet dataSet)
        throws ObjectStoreException {
        String header = ((DNASequence) bioJavaSequence).getOriginalHeader();
        String mrnaIdentifier = getIdentifier(bioJavaSequence);
        ObjectStore os = getIntegrationWriter().getObjectStore();
        Model model = os.getModel();
        if (model.hasClassDescriptor(model.getPackageName() + ".CDS")) {
            Class<? extends FastPathObject> cdsCls =
                model.getClassDescriptorByName("CDS").getType();
            if (!DynamicUtil.isInstance(bioEntity, cdsCls)) {
                throw new RuntimeException("the InterMineObject passed to "
                        + "MedicagoCDSFastaLoaderTask.extraProcessing() is not a "
                        + "CDS: " + bioEntity);
            }
            InterMineObject mrna = getMRNA(mrnaIdentifier, organism, model);
            if (mrna != null) {
                bioEntity.setFieldValue("transcript", mrna);
            }
            Location loc = getLocationFromHeader(header, (SequenceFeature) bioEntity,
                    organism);
            getDirectDataLoader().store(loc);
        } else {
            throw new RuntimeException("Trying to load CDS sequence but CDS does not exist in the"
                    + " data model");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    protected String getIdentifier(Sequence bioJavaSequence) {
        String header = bioJavaSequence.getAccession().getID();
        final String regexp = "^(\\S+)\\s+\\|.+";
        Pattern p = Pattern.compile(regexp);
        Matcher m = p.matcher(header);
        if (m.matches()) {
          String mrnaIdentifier = m.group(1);
          return mrnaIdentifier;
        }
        throw new RuntimeException("header doesn't match pattern \"" + regexp + "\": " + header);
    }
}
