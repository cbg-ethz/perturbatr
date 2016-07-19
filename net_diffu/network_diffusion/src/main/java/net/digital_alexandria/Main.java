package net.digital_alexandria;


import net.digital_alexandria.mvc.view.GraphForm;

public class Main
{

    public static void main(String[] args)
    {
        if (args.length != 0 && args[0].equals("--gui"))
            GraphForm.run(args);
        else
            new Controller(args).run();
    }

}
