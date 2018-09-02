classdef Season
    enumeration
        Winter, Spring, Summer, Fall
    end
    methods(Static)
        function season = getSeason(month)
            if month == 1 || month ==2 || month == 12
                season = 'Winter';
            elseif month == 3 || month ==4 || month ==5
                season = 'Spring';
            elseif month == 6 || month ==7 || month ==8
                season = 'Summer';
            elseif month == 9 || month ==10 || month ==11
                season = 'Fall';
            end
        end
        function marker = getSeasonMarker(season)
            % Winter
            if season == 1
                marker = 'o';
            % Spring                
            elseif season == 2
                marker = '+';
            % Summer   
            elseif season == 3
                marker = '*';
            % Fall
            elseif season == 4
                marker = '.';
            end
        end   
%         function season = getSeason(month)
%             if month == 1 || month ==2 || month == 12
%                 season = 1;
%             elseif month == 3 || month ==4 || month ==5
%                 season = 2;
%             elseif month == 6 || month ==7 || month ==8
%                 season = 3;
%             elseif month == 9 || month ==10 || month ==11
%                 season = 4;
%             end
%         end
    end
end